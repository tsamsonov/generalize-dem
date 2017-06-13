import os
import sys
import arcpy
import multiprocessing


def worker(oid):
    arcpy.AddMessage(oid)
    return


def execute(*oids):
    # Set multiprocessing exe in case we're running as an embedded process, i.e ArcGIS
    # get_install_path() uses a registry query to figure out 64bit python exe if available

    arcpy.AddMessage('Setting')
    multiprocessing.set_executable(os.path.join(get_install_path(), 'pythonw.exe'))

    arcpy.AddMessage('Pooling')
    pool = multiprocessing.Pool()

    arcpy.AddMessage('Working')
    pool.map(worker, oids)

    return

def get_install_path():
    ''' Return 64bit python install path from registry (if installed and registered),
        otherwise fall back to current 32bit process install path.
    '''
    if sys.maxsize > 2**32: return sys.exec_prefix #We're running in a 64bit process

    # We're 32 bit so see if there's a 64bit install
    path = r'SOFTWARE\Python\PythonCore\2.7'

    from _winreg import OpenKey, QueryValue
    from _winreg import HKEY_LOCAL_MACHINE, KEY_READ, KEY_WOW64_64KEY

    try:
        with OpenKey(HKEY_LOCAL_MACHINE, path, 0, KEY_READ | KEY_WOW64_64KEY) as key:
            return QueryValue(key, "InstallPath").strip(os.sep) # We have a 64bit install, so return that.
    except: return sys.exec_prefix # No 64bit, so return 32bit path


if __name__=='__main__':
    import Test
    oids = arcpy.GetParameterAsText(0)

    Test.execute(oids)
