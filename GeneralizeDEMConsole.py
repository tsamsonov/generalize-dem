import GeneralizeDEM
import multiprocessing
import arcpy
import time

def worker(oid):
    arcpy.AddMessage('Yeah')
    time.sleep(5)
    return oid

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

    # arcpy.AddMessage('Trying to make multiprocessing')
    # arcpy.AddMessage('Creating Pool')
    # pool = multiprocessing.Pool(processes=3)
    # jobs = []
    # oids = [1, 2, 3, 4, 5]
    # arcpy.AddMessage(multiprocessing.cpu_count())

    # PROCESS
    # try:
    #     for oid in oids:
    #         p = multiprocessing.Process(target = worker, args = (oid,))
    #         jobs.append(p)
    #         p.start()
    #     for job in jobs:
    #         job.join()
    # except:
    #     input("Press Enter to continue...")

    # POOL ASYNC
    # try:
    #     for oid in oids:
    #         pool.apply_async(worker, (oid,))
    #     pool.close()
    #     pool.join()
    # except:
    #     input("Press Enter to continue...")

    # POOL MAP
    # try:
    #     pool.map(worker, oids)
    #     pool.close()
    #     pool.join()
    # except:
    #     input("Press Enter to continue...")


    input("Press Enter to continue...")