import ast
import xml.etree.ElementTree as ET
import base64
from datetime import datetime, timedelta
import pytz
import struct
from unidecode import unidecode

def create_jcamp_library(df):
    jcamp_library = []
    
    for index, row in df.iterrows():
        if type(row['MS_Peaks']) == str:
            row['MS_Peaks'] = ast.literal_eval(row['MS_Peaks'].strip())
        jcamp_entry = f"""##TITLE={row['Chemical_Name']}
##JCAMPDX=Revision 5.00
##DATA TYPE=MASS SPECTRUM
##ORIGIN=Combined CEF Results
##OWNER=
##CAS REGISTRY NO={row['CAS_Number']}
##$RETENTION INDEX={row['RI']}
##RETENTION TIME={row['RT']}
##MOLECULAR FORMULA={row['Formula']}
##XUNITS=M/Z
##YUNITS=RELATIVE INTENSITY
##XFACTOR=1
##YFACTOR=1
##FIRSTX=0
##LASTX={max(peak[0] for peak in row['MS_Peaks'])}
##FIRSTY=0
##MAXX={max(peak[0] for peak in row['MS_Peaks'])}
##MINX=0
##MAXY=999
##MINY=0
##MW={row['RT']}
##NPOINTS={len(row['MS_Peaks'])}
##PEAK TABLE=(XY..XY)
"""     

        for mz, intensity in row['MS_Peaks']:
            jcamp_entry += f"{int(mz)} {int(intensity)}\n"
        
        jcamp_entry += "##END="
        jcamp_library.append(jcamp_entry)
    
    return jcamp_library

def write_jcamp_library(df, filename):
    jcamp_library = create_jcamp_library(df)
    with open(filename, 'w') as f:
        for entry in jcamp_library:
            f.write(entry + '\n\n')

def get_custom_datetime_string():
    # Get current time with timezone UTC
    now = datetime.now(pytz.timezone('UTC'))
    
    # Convert it to a specific timezone (in this case, UTC+12:00)
    target_timezone = pytz.timezone('Pacific/Auckland')  # +12:00 timezone
    now_in_tz = now.astimezone(target_timezone)

    # Format the date-time string with 7 digits of precision for microseconds
    formatted_time = now_in_tz.strftime('%Y-%m-%dT%H:%M:%S.%f')
    
    # Manually format microseconds to 7 digits instead of the default 6
    microseconds_str = str(now_in_tz.microsecond).ljust(7, '0')
    
    # Add the timezone offset in '+HH:MM' format
    timezone_offset = now_in_tz.strftime('%z')
    timezone_offset_formatted = f'{timezone_offset[:3]}:{timezone_offset[3:]}'

    # Combine all parts into the final format
    custom_datetime_string = f"{formatted_time[:-6]}{microseconds_str}{timezone_offset_formatted}"
    
    return custom_datetime_string

def create_mslibrary_xml(df, filename):
    root = ET.Element("LibraryDataSet", SchemaVersion="2", xmlns="Quantitation.LibraryDatabase")
    edit_datetime = get_custom_datetime_string()
    library = ET.SubElement(root, "Library")
    ET.SubElement(library, "LibraryID").text = "1"
    ET.SubElement(library, "CreationDateTime").text = edit_datetime

    
    for index, row in df.iterrows():
        compound = ET.SubElement(root, "Compound")
        ET.SubElement(compound, "LibraryID").text = "1"
        ET.SubElement(compound, "CompoundID").text = str(index + 1)
        ET.SubElement(compound, "AlternateNames").text = ""  # You may want to add this to your DataFrame
        ET.SubElement(compound, "CASNumber").text = row['CAS_Number']
        ET.SubElement(compound, "CompoundName").text = unidecode(row['Chemical_Name'])
        ET.SubElement(compound, "Description").text = row['File']  # You may want to add this to your DataFrame
        ET.SubElement(compound, "Formula").text = row['Formula']
        ET.SubElement(compound, "LastEditDateTime").text = edit_datetime
        ET.SubElement(compound, "RetentionTimeRTL").text = str(row['RT'])  # Using RT as MW, adjust if needed
        ET.SubElement(compound, "RetentionIndex").text = str(row['RI'])
        ET.SubElement(compound, "UserDefined").text = str(row['RI Ref'])  # You may want to add this to your DataFrame

        spectrum = ET.SubElement(root, "Spectrum")
        ET.SubElement(spectrum, "LibraryID").text = "1"
        ET.SubElement(spectrum, "CompoundID").text = str(index + 1)
        ET.SubElement(spectrum, "SpectrumID").text = "0"
        if type(row['MS_Peaks']) is str:
            ms_peaks = ast.literal_eval(row['MS_Peaks'].strip())
        else:
            ms_peaks = row['MS_Peaks']
        mz_values = [peak[0] for peak in ms_peaks]
        abundance_values = [peak[1] for peak in ms_peaks]

        ET.SubElement(spectrum, "AbundanceValues").text = base64.b64encode(struct.pack('<%sd' % len(abundance_values), *abundance_values)).decode()
        ET.SubElement(spectrum, "IonPolarity").text = "Positive"
        ET.SubElement(spectrum, "LastEditDateTime").text = edit_datetime
        ET.SubElement(spectrum, "MzValues").text = base64.b64encode(struct.pack('<%sd' % len(mz_values), *mz_values)).decode()
        ET.SubElement(spectrum, "Origin").text = "Combined CEF Results"
        ET.SubElement(spectrum, "Owner").text = ""
        ET.SubElement(spectrum, "ScanType").text = "Scan"

    tree = ET.ElementTree(root)
    # Convert the ElementTree to a string
    xml_string = ET.tostring(root, encoding="utf-8", xml_declaration=True).decode("utf-8")
    
    # Prettify the XML string
    from xml.dom import minidom
    pretty_xml_string = minidom.parseString(xml_string).toprettyxml(indent="  ")
    
    # Return the prettified XML string
    return pretty_xml_string
    # print(filename)
    # tree.write(filename, encoding="utf-8", xml_declaration=True)
