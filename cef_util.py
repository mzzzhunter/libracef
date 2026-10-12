import pandas as pd
import xml.etree.ElementTree as ET
import os
import re
import numpy as np
from scipy.spatial.distance import cosine

def parse_cef_file(file_path):
    tree = ET.parse(file_path)
    root = tree.getroot()

    compounds = []
    for compound in root.findall('.//Compound'):
        compound_data = {}
        
        # Find the chemical name and formula in the attributes under Molecule node under Results node
        results_node = compound.find('.//Results')
        if results_node is not None:
            molecule_node = results_node.find('.//Molecule')
            if molecule_node is not None:
                compound_data['Chemical_Name'] = molecule_node.get('name')
                compound_data['Formula'] = molecule_node.get('formula')
                
                # Get CAS number from Accession node under Database node under Molecule node
                database_node = molecule_node.find('.//Database')
                if database_node is not None:
                    accession_node = database_node.find('.//Accession')
                    if accession_node is not None:
                        compound_data['CAS_Number'] = accession_node.get('id')
        
        # Find RT and RI in attributes in location node under Compound node
        location_node = compound.find('.//Location')
        if location_node is not None:
            compound_data['RT'] = float(location_node.get('rt'))
            compound_data['RI'] = float(location_node.get('ri', 0))  # Default to 0 if 'ri' is not present
            compound_data['MaxArea'] = float(location_node.get('a', 0))  # Default to 0 if 'area' is not present

        # Capture mass spectrum from MSPeaks node under Spectrum node
        spectrum_node = compound.find('.//Spectrum')
        if spectrum_node is not None:
            ms_peaks = spectrum_node.find('.//MSPeaks')
            if ms_peaks is not None:
                peaks = []
                for peak in ms_peaks.findall('.//p'):
                    peaks.append((round(float(peak.get('x'))), float(peak.get('y'))))
                compound_data['MS_Peaks'] = peaks

        
        compound_data['File'] = os.path.basename(file_path)  # Add filename to track source
        compounds.append(compound_data)

    return compounds

def blank_compound_mask(data, blank_compounds):
    """Match CAS when both records have one, otherwise match nonempty names."""
    def text(value):
        return '' if pd.isna(value) else str(value).strip()

    def name(value):
        return re.sub(r'\s+peak\s+(?:[2-9]|[1-9]\d+)$', '', text(value).casefold())

    blank_cas = set()
    blank_names = set()
    names_without_cas = set()
    for compound in blank_compounds:
        cas = text(compound.get('CAS_Number'))
        chemical_name = name(compound.get('Chemical_Name'))
        if cas:
            blank_cas.add(cas)
        if chemical_name:
            blank_names.add(chemical_name)
            if not cas:
                names_without_cas.add(chemical_name)

    matches = []
    for _, compound in data.iterrows():
        cas = text(compound.get('CAS_Number'))
        chemical_name = name(compound.get('Chemical_Name'))
        matches.append(bool(cas and cas in blank_cas) or
                       bool(chemical_name) and chemical_name in
                       (names_without_cas if cas else blank_names))
    return pd.Series(matches, index=data.index, dtype=bool)


def combine_cef_results(directory, rt_tolerance=0.1, group_similarity_threshold=0.9, exclude_elements=['Si']):
    all_compounds = []
    for filename in os.listdir(directory):
        if filename.endswith('.cef'):
            file_path = os.path.join(directory, filename)
            all_compounds.extend(parse_cef_file(file_path))
    
    # Sort compounds by Chemical_Name and RT
    sorted_list = sorted(all_compounds, key=lambda x: (x['Chemical_Name'], x['RT']))

    # Merge adjacent same-name observations within the RT tolerance, retaining
    # the highest-area observation independently of the previous input record.
    unique_compounds = []
    previous = None
    for current in sorted_list:
        if (previous is None
                or current['Chemical_Name'] != previous['Chemical_Name']
                or abs(current['RT'] - previous['RT']) >= rt_tolerance):
            unique_compounds.append(current.copy())
        elif current['MaxArea'] > unique_compounds[-1]['MaxArea']:
            retained = unique_compounds[-1]
            for field in ('MaxArea', 'MS_Peaks', 'File', 'RT', 'RI'):
                retained[field] = current[field]
        previous = current

    # Step 3: Rename duplicates with "peak" suffix
    # Create a new list to hold renamed compounds
    renamed_dicts = []
    peak_counts = {}  # To keep track of peak counts for renaming

    for peak in unique_compounds:
        if peak['Chemical_Name'] not in peak_counts:
            peak_counts[peak['Chemical_Name']] = 1
        else:
            peak_counts[peak['Chemical_Name']] += 1

        # If there are multiple peaks with the same name, rename them
        if peak_counts[peak['Chemical_Name']] > 1:
            renamed_peak = peak.copy()
            renamed_peak['Chemical_Name'] = f"{peak['Chemical_Name']} peak {peak_counts[peak['Chemical_Name']]}"
            renamed_dicts.append(renamed_peak)
        else:
            renamed_dicts.append(peak)

    # Sort unique_compounds by RT
    renamed_dicts.sort(key=lambda x: x['RT'])

    column_order = ['Chemical_Name', 'Formula', 'RT', 'RI', 'RI Ref', 'CAS_Number', 'group', 'Similarity_to_Previous', 'Similarity_to_Next', 'MS_Peaks', 'File','MaxArea']
    df = pd.DataFrame(renamed_dicts)
    if df.empty:
        return pd.DataFrame(columns=column_order)
    
    # Ignore blank entries from the comma-separated exclusion input.
    for element in exclude_elements:
        element = element.strip()
        if element:
            df = df[~df['Formula'].str.contains(element, na=False)]

    # Group assignment below uses index label 0 for the first retained row.
    df = df.reset_index(drop=True)
    if df.empty:
        return pd.DataFrame(columns=column_order)
    
    # Calculate similarity to previous and next spectrum
    df['Similarity_to_Previous'] = df.apply(lambda row: calculate_similarity(row, df.iloc[df.index.get_loc(row.name) - 1] if df.index.get_loc(row.name) > 0 else None), axis=1)
    df['Similarity_to_Next'] = df.apply(lambda row: calculate_similarity(row, df.iloc[df.index.get_loc(row.name) + 1] if df.index.get_loc(row.name) < len(df) - 1 else None), axis=1)
    # Replace None values with 0 in Similarity_to_Previous and Similarity_to_Next columns
    df['Similarity_to_Previous'] = df['Similarity_to_Previous'].fillna(0)
    df['Similarity_to_Next'] = df['Similarity_to_Next'].fillna(0)
    # Add a 'group' column and assign group numbers
    df['group'] = 0
    group_id = 1

    # Assign the first row to the first group
    df.at[0, 'group'] = group_id

    for i in range(1, len(df)):
    # Since next_similarity of previous row == prev_similarity of current row,
    # we check if this shared similarity exceeds the threshold
        shared_similarity = df.iloc[i]['Similarity_to_Previous']  # or df.at[i - 1, 'next_similarity']

        if shared_similarity > group_similarity_threshold:
            # Assign the current row to the same group as the previous row
            df.iloc[i, df.columns.get_loc('group')] = group_id
        else:
            # Start a new group
            group_id += 1
            df.iloc[i, df.columns.get_loc('group')] = group_id

    df['RI Ref'] = ''
    # Reorder columns
    df = df[column_order]
    
    return df

def calculate_similarity(row1, row2):
    if row2 is None:
        return None
    
    peaks1 = row1['MS_Peaks']
    peaks2 = row2['MS_Peaks']
    
    # Create vectors for comparison
    max_mz = max(max(p[0] for p in peaks1), max(p[0] for p in peaks2))
    vector1 = np.zeros(int(max_mz) + 1)
    vector2 = np.zeros(int(max_mz) + 1)
    
    for mz, intensity in peaks1:
        vector1[int(mz)] = intensity
    for mz, intensity in peaks2:
        vector2[int(mz)] = intensity
    
    # Calculate cosine similarity
    similarity = 1 - cosine(vector1, vector2)
    return similarity

# # Usage
# cef_directory = 'path/to/your/cef/files/'
# combined_results_df = combine_cef_results(cef_directory)

# # Now you can work with the combined results in a pandas DataFrame
# print(combined_results_df.head())
# print(f"Total compounds found: {len(combined_results_df)}")















