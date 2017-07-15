import sys
import time
import arcpy
import traceback
import GeneralizeDEM

if __name__ == '__main__':

    # SET PARAMETERS HERE
    # --------------------------------------------------------------------
    demdataset = 'X:/Work/DEMGEN/chukotka_Clip'
    marine = 'X:/Work/Scripts & Tools/MY/DEMGEN/DEMGENEW.gdb/ne_10m_ocean_P'
    output = 'X:/Work/DEMGEN/DEMGENEW.gdb/chukotka_gen'
    outputcellsize = 2000
    minacc1 = 40
    minlen1 = 10
    minacc2 = 20
    minlen2 = 5
    is_widen = 'true'
    widentype = 'Min/Max'
    widendist = 4000
    filtersize = 5
    is_smooth = 'true'
    is_parallel = 'true'
    tilesize = 1100
    # --------------------------------------------------------------------

    print('> Initializing GeneralizeDEM script...')
    print('')

    start = int(time.time())
    try:
        GeneralizeDEM.execute(demdataset, marine, output, outputcellsize,
                              minacc1, minlen1, minacc2, minlen2,
                              is_widen, widentype, widendist, filtersize,
                              is_smooth, tilesize, 0, is_parallel)
    except Exception:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        pymsg = "Traceback Info:\n" + tbinfo + "\nError Info:\n    " + \
                str(sys.exc_type) + ": " + str(sys.exc_value) + "\n"
        arcpy.AddError(pymsg)
        print("Processing failed")

    finish = int(time.time())
    seconds = finish - start
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)

    print ''
    print "> Finished in %02d h %02d m %02d s" % (h, m, s)
    print ''

    input("Press Enter to continue...")