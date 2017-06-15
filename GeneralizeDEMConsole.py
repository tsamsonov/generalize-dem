import GeneralizeDEM

if __name__ == '__main__':

    demdataset = 'X:/Work/Scripts & Tools/MY/DEMGEN/mistral'
    marine = None
    output = 'X:/Work/Scripts & Tools/MY/DEMGEN/DEMGENEW.gdb/parallel'
    outputcellsize = 1000
    minacc1 = 40
    minlen1 = 20
    minacc2 = 20
    minlen2 = 10
    is_widen = 'true'
    widentype = 0
    widendist = 0
    filtersize = 0
    is_smooth = 'true'
    is_parallel = 'true'
    tilesize = 512

    try:
        GeneralizeDEM.execute(demdataset, marine, output, outputcellsize,
                              minacc1, minlen1, minacc2, minlen2,
                              is_widen, widentype, widendist, filtersize,
                              is_smooth, is_parallel, tilesize)
    except Exception:
        print("Quit with error...")
        input("Press Enter to continue...")


    input("Press Enter to continue...")