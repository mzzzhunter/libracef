import ast
import subprocess
# NIST search
def format_ms_peaks_for_nist(peaks):
    if type(peaks) == str:
        peaks = ast.literal_eval(peaks)
    return ''.join([f"{int(mz)}\t{int(intensity)}\n" for mz, intensity in peaks])
    
def create_spectrum_file(df_row, filename):
    lines = []
    spec_peaks = format_ms_peaks_for_nist(df_row['MS_Peaks'])
    lines.append(f"name:RT - {df_row['RT']}\n")
    lines.append(f"num:{len(df_row['MS_Peaks'])}\n")
    lines.append(spec_peaks)
    
    with open(filename, 'w') as f:
        f.writelines(lines)
    return filename

def get_autoimp_path(nist_path="c:\\NIST14\\MSSEARCH\\"):
    with open(nist_path+"AUTOIMP.MSD") as f:
        path = f.readline()
    return path.strip()

def update_filespec(filespec_path, data_path):
    with open(filespec_path, 'w') as f:
        f.write(data_path + '\n10 724')


def search_nist_for_spectrum(df_row, spec_data_path='c:\\NIST14\\MSSEARCH\\datafile.txt', nist_path=r"C:\NIST14\MSSEARCH\nistms$.exe"):
    spec_data_file = create_spectrum_file(df_row, spec_data_path)
    filespec_path = get_autoimp_path()
    update_filespec(filespec_path, spec_data_file)
    cmd = [
        nist_path,
        "/INSTRUMENT"
    ]
    try:
        # Execute the NIST MS Search command
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running NIST MS Search: {e}")