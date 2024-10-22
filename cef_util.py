import pandas as pd
import xml.etree.ElementTree as ET
import os
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

def combine_cef_results(directory, rt_tolerance=0.1, group_similarity_threshold=0.9, exclude_elements=['Si']):
    all_compounds = []
    for filename in os.listdir(directory):
        if filename.endswith('.cef'):
            file_path = os.path.join(directory, filename)
            all_compounds.extend(parse_cef_file(file_path))
    
    # Sort compounds by Chemical_Name and RT
    sorted_list = sorted(all_compounds, key=lambda x: (x['Chemical_Name'], x['RT']))

    # Step 2: Prepare the structures to hold the unique peaks and count the duplicates
    unique_compounds = []
    peak_counts = {}  # To keep track of peak counts for renaming

    for i, current in enumerate(sorted_list):
        if i == 0:
            # Always add the first peak
            unique_compounds.append(current)
            
        else:
            previous = sorted_list[i - 1]
            
            # Check if the current peak is not a duplicate based on name and RT difference
            if current['Chemical_Name'] != previous['Chemical_Name'] or abs(current['RT'] - previous['RT']) >= 0.1:
                # It's not a duplicate or it has a significant RT difference
                unique_compounds.append(current)

    # Step 3: Rename duplicates with more than 0.1 RT difference
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

    # Rename duplcates as Peak 1, Peak 2, etc.
    
    df = pd.DataFrame(renamed_dicts)
    
    # Remove compounds with formula containing "Si"
    for element in exclude_elements:
        df = df[~df['Formula'].str.contains(element, na=False)]
    
    
    # Calculate similarity to previous and next spectrum
    df['Similarity_to_Previous'] = df.apply(lambda row: calculate_similarity(row, df.iloc[df.index.get_loc(row.name) - 1] if df.index.get_loc(row.name) > 0 else None), axis=1)
    df['Similarity_to_Next'] = df.apply(lambda row: calculate_similarity(row, df.iloc[df.index.get_loc(row.name) + 1] if df.index.get_loc(row.name) < len(df) - 1 else None), axis=1)
    # Replace None values with 0 in Similarity_to_Previous and Similarity_to_Next columns
    df['Similarity_to_Previous'] = df['Similarity_to_Previous'].fillna(0)
    df['Similarity_to_Next'] = df['Similarity_to_Next'].fillna(0)
    # Add a 'group' column and assign group numbers
    df['group'] = 0
    current_group = 1
    
    for i in range(len(df)):
        if i == 0:
            df.iloc[i, df.columns.get_loc('group')] = current_group
        else:
            prev_similarity = df.iloc[i]['Similarity_to_Previous']
            next_similarity = df.iloc[i]['Similarity_to_Next']
            
            if (prev_similarity is not None and prev_similarity > group_similarity_threshold) or (next_similarity is not None and next_similarity > group_similarity_threshold):
                if i > 0 and df.iloc[i-1]['group'] != 0:
                    df.iloc[i, df.columns.get_loc('group')] = df.iloc[i-1]['group']
                else:
                    df.iloc[i, df.columns.get_loc('group')] = current_group
            else:
                current_group += 1
                df.iloc[i, df.columns.get_loc('group')] = current_group
    
    # Clean up any remaining ungrouped rows
    for i in range(len(df)):
        if df.iloc[i]['group'] == 0:
            current_group += 1
            df.iloc[i, df.columns.get_loc('group')] = current_group
    
    # Ensure 'group' is included in the column order
    # column_order.insert(-1, 'group')  # Insert 'group' before the last column ('File')
    df['RI Ref'] = ''
    # Reorder columns
    column_order = ['Chemical_Name', 'Formula', 'RT', 'RI', 'RI Ref', 'CAS_Number', 'group', 'Similarity_to_Previous', 'Similarity_to_Next', 'MS_Peaks', 'File']
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















