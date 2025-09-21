import subprocess 
import os
import argparse
 
if __name__ == "__main__":
    chemin_dossier = "dsirp/standard/"
    listInstance = []
    liste = ["25-10-3",
            "100-20-4",
            "125-20-1",
            "125-20-2",
            "125-20-3",
            "125-20-4",
            "125-20-5",
            "150-20-1",
            "150-20-2",
            "150-20-3",
            "150-20-4",
            "150-20-5",
            "200-20-1",
            "200-20-2",
            "200-20-3",
            "200-20-4",
            "200-20-5",
            "200-10-1",
            "200-10-2",
            "200-10-3",
            "200-10-4",
            "200-10-5"]
    # for nbClients in [5,10,15,25,50,75,100,125,150,200]:
    #     for p in [5,10,20]:
    #         for i in [1,2,3,4,5]:
    for s in liste:
        # listInstance.append("dirp-" + str(nbClients) + "-" + str(p) + "-" + str(i) + ".dat")
        listInstance.append("dirp-" + s + ".dat" )

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
        chemin_element = os.path.join(chemin_dossier, listInstance[expe_id%150])
        deterministic = 0
        endDay = 0

        if (True):
            nbScenario = 1
            maxIter = 80000
            maxIterNonProd = 80000
            maxTime = 1800
            deterministic = 1
            endDay = 1

        if os.path.isfile(chemin_element):
            result = subprocess.Popen(['./irp', chemin_element, '-veh', str(nbVehicle), '-seed', str(seed), '-scenarios', str(nbScenario), '-iter', str(maxIter), "-iterNonProd", str(maxIterNonProd), '-time', str(maxTime), '-traces', '0', '-cores', str(nbCores), '-endDay', str(endDay), '-deterministic', str(deterministic)], stdout=subprocess.PIPE, text=True, stderr=subprocess.PIPE)
            for res in result.stdout:
                print(res)
            for err in result.stderr:
                print(err)