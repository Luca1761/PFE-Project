import subprocess 
import os
import argparse
 
if __name__ == "__main__":
    chemin_dossier = "dsirp/standard/"
    listInstance = []
    for nbClients in [5,10,15,25,50,75,100,125,150,200]:
        for p in [5,10,20]:
            for i in [1,2,3,4,5]:
                listInstance.append("dirp-" + str(nbClients) + "-" + str(p) + "-" + str(i) + ".dat")

    launchAll = False

    if launchAll:
        average = 0
        errorList = []
        for nom_element in listInstance:
            chemin_element = os.path.join(chemin_dossier, nom_element)
            if os.path.isfile(chemin_element):
                print(chemin_element)
                result = subprocess.Popen(['./irp.exe', chemin_element, '-veh', '1',  '-seed',  '55',  '-scenarios', '1', '-iter', '10000', "-iterNonProd", '10000', '-time', '3', '-traces', '0', '-cores', '4'], stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE)
                for res in result.stdout:
                    average += float(res)
                    print(res)
                for err in result.stderr:
                    errorList.append([chemin_element, err])
                    print(err)
        print(average/len(listInstance) )

    else:
        parser = argparse.ArgumentParser(description = 'Launch experiences')
        parser.add_argument('--expe_id', type = int, default = 0)
        args = parser.parse_args()

        expe_id = args.expe_id
        nbScenario = 10
        nbVehicle = 1
        seed = 42
        nbCores = 10
        maxIter = 20000
        maxIterNonProd = 20000
        maxTime = 1200
        chemin_element = os.path.join(chemin_dossier, listInstance[expe_id])

        if os.path.isfile(chemin_element):
            result = subprocess.Popen(['./irp', chemin_element, '-veh', str(nbVehicle), '-seed', str(seed), '-scenarios', str(nbScenario), '-iter', str(maxIter), "-iterNonProd", str(maxIterNonProd), '-time', str(maxTime), '-traces', '0', '-cores', str(nbCores), '-endDay', '0', '-deterministic', '0'], stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE)
            for res in result.stdout:
                print(res)
            for err in result.stderr:
                print(err)