import subprocess 
import os
import argparse
 
if __name__ == "__main__":
    listInstance = []
    for nbClients in [5,10,15,25,50,75,100,125,150,200]:
        for p in [5,10,20]:
            for i in [1,2,3,4,5]:
                listInstance.append("dirp-" + str(nbClients) + "-" + str(p) + "-" + str(i) + ".dat")

    parser = argparse.ArgumentParser(description = 'Launch experiences')
    parser.add_argument('--expe_id', type = int, default = 0)
    parser.add_argument('--instance_type', type = str, default = "seasonal")
    args = parser.parse_args()

    expe_id = args.expe_id
    instance_type = args.instance_type
    folder_path = "instances/" + instance_type + "/"
    file_path = os.path.join(folder_path, listInstance[expe_id])

    nbScenario = 10
    nbVehicle = 1
    seed = 42
    maxIter = 20000
    maxIterNonProd = 20000
    maxTime = 1200
    deterministic = 0
    endDayInventories = 0

    if os.path.isfile(file_path):
        result = subprocess.Popen(['./irp', file_path, '-veh', str(nbVehicle), 
                                   '-seed', str(seed), '-scenarios', str(nbScenario), 
                                   '-iter', str(maxIter), "-iterNonProd", str(maxIterNonProd), 
                                   '-time', str(maxTime), '-traces', '0', '-cores', str(nbScenario), 
                                   '-endDay', str(endDayInventories), '-deterministic', 
                                   str(deterministic)], stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE)
        for res in result.stdout:
            print(res)
        for err in result.stderr:
            print(err)