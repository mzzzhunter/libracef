import ast
import subprocess
import pathlib
from PyQt5.QtWidgets import QMessageBox

# NIST search
def format_ms_peaks_for_nist(peaks):
    if type(peaks) == str:
        peaks = ast.literal_eval(peaks)
    return ''.join([f"{int(mz)}\t{int(intensity)}\n" for mz, intensity in peaks])
    
def create_spectrum_file(df_row, filename):
    lines = []
    peaks = df_row['MS_Peaks']
    if isinstance(peaks, str):
        peaks = ast.literal_eval(peaks.strip())
    spec_peaks = format_ms_peaks_for_nist(peaks)
    lines.append(f"name:RT - {df_row['RT']}\n")
    lines.append(f"num:{len(peaks)}\n")
    lines.append(spec_peaks)
    
    with open(filename, 'w') as f:
        f.writelines(lines)
    return filename

def get_autoimp_path(nist_path):
    nist_directory = pathlib.Path(nist_path)
    autoimp_path = nist_directory / "AUTOIMP.MSD"

    try:
        with autoimp_path.open() as f:
            path = f.readline()
    except FileNotFoundError:
        filespec_path = nist_directory / "FILESPEC.FIL"
        autoimp_path.write_text(str(filespec_path))
        return str(filespec_path)

    return path.strip()

def check_autoimp_writable(nist_path):
    autoimp_path = pathlib.Path(nist_path) / "AUTOIMP.MSD"
    if not autoimp_path.exists():
        get_autoimp_path(nist_path)

    with autoimp_path.open("r+"):
        pass


def update_filespec(filespec_path, data_path):
    with open(filespec_path, 'w') as f:
        f.write(data_path + '\n10 724')


def show_nist_error(parent, message):
    QMessageBox.critical(parent, "NIST MS Search Error", message)


def search_nist_for_spectrum(df_row, spec_data_path, nist_path, parent=None):
    try:
        spec_data_file = create_spectrum_file(df_row, spec_data_path)
    except (OSError, KeyError, TypeError, ValueError) as e:
        show_nist_error(parent, f"Could not create the NIST spectrum file:\n{e}")
        return

    try:
        filespec_path = get_autoimp_path(nist_path)
    except OSError as e:
        show_nist_error(parent, f"Could not read or create NIST AUTOIMP.MSD:\n{e}")
        return

    try:
        update_filespec(filespec_path, spec_data_file)
    except OSError as e:
        show_nist_error(parent, f"Could not update the NIST filespec:\n{e}")
        return

    cmd = [
        nist_path+"\\nistms$.exe",
        "/INSTRUMENT"
    ]
    try:
        # Execute the NIST MS Search command
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        details = e.stderr.strip() if e.stderr else str(e)
        show_nist_error(parent, f"NIST MS Search failed to run:\n{details}")
    except OSError as e:
        show_nist_error(parent, f"Could not launch NIST MS Search:\n{e}")
