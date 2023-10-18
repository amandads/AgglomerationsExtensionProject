from scipy.stats import mannwhitneyu
from cliffs_delta import cliffs_delta
import subprocess
import pandas as pd


metrics = ["wmc","maxNest","cbo","lcom*","noc","dit"]
listaGeneral = ["heterogeneous","aggFE","aggLM","uniqueLC","uniqueRB","uniqueFE","uniqueLM","wbs"]
listaDeep = ["FeLm","LcFe","LcLm","LcLmFe","LcRb","LcRbFe","LcRbLm","LcRbLmFe","RbFe","RbLm","RbLmFe"]
listaDeep2 = ["FeLm","LcFe","LcLm","LcLmFe","LcRb","LcRbLm","LcRbLmFe","RbFe","RbLm","RbLmFe"]
myPathGEQC = ""
myPathDEQC = ""
myPathGEGH = ""
myPathDEGH = ""
    

def calculaMan(path,metric,position,meuArqX,meuArqY,minhaSaida):

    if metric=="lcom*":
        opt = metric.strip("*")
    else:
        opt = metric
        
    myStringDir = "mkdir \""+path+opt
    #print(myStringDir)
    makingDir = subprocess.run(myStringDir,shell=True)

    arquivoX = open(path+meuArqX+".txt","r")
    
    arquivoY = open(path+meuArqY+".txt","r")
    headerX = arquivoX.readline()
    headerY = arquivoY.readline()
    saida = open(path+opt+"/"+minhaSaida+".txt","w")


    meuX = [0]
    meuY = [0]

    for lines in arquivoX:
        if ";NaN;" not in lines:
            myString = lines.split(";")
            meuX.append(float(myString[position])+0.000000000000001)


    for lines in arquivoY:
        if ";NaN;" not in lines:        
            myString = lines.split(";")
            meuY.append(float(myString[position])+0.000000000000001)


    #print(meuX)
    
    U1, p = mannwhitneyu(meuX,meuY,method="auto")

    saida.write(meuArqX+";"+meuArqY+":     "+str(U1)+"    "+str(p))


    arquivoX.close()
    arquivoY.close()
    saida.close()



def funcaoIteradora(path, lista):

    cont = 8
    for items in metrics:
        #print("METRICS: "+items+"\n")
        for elements in lista:
            for nomes in lista:
                if elements!=nomes:
                    #print(elements+"    "+nomes+"\n")
                    calculaMan(path,items,cont,elements,nomes,elements+"-"+nomes)

        #print("\n\n")
        cont = cont + 1
    


funcaoIteradora(myPathGEQC,listaGeneral)
funcaoIteradora(myPathDEQC,listaDeep)
funcaoIteradora(myPathGEGH,listaGeneral)
funcaoIteradora(myPathDEGH,listaDeep2)

print("Finished U")

def myBonferroni(numeroDivisor,meuPValue,arqE,arqS,name):

    minhaComp = meuPValue / numeroDivisor
    #print(minhaComp)


    entrada = open(arqE,"r")
    saida = open(arqS,"a")

    myLine = entrada.readline()
    myString = myLine.split("    ")
    pvalue = float(myString[2])
    #print(pvalue)

    if(pvalue<=minhaComp):
        #print(name+";SUCCESS;"+str(pvalue)+";"+str(minhaComp)+"\n")
        saida.write(name+";SUCCESS;"+str(pvalue)+";"+str(minhaComp)+"\n")
    else:
        #print(name+";FAILED;"+str(pvalue)+";"+str(minhaComp)+"\n")
        saida.write(name+";FAILED;"+str(pvalue)+";"+str(minhaComp)+"\n")

    entrada.close()
    saida.close()

def secondIterator(number,pvalue,path,metric,lista):

    for elements in metric:
        if elements == "lcom*":
            opt = elements.strip("*")
        else:
            opt = elements
        saida = path+opt+"/bonferroni.txt"
        cleanFile = open(saida,"w")
        cleanFile.close()
        
        for items in lista:

            for keys in lista:

                if items!=keys:
                    pathXY = path+opt+"/"+items+"-"+keys+".txt"
                    myBonferroni(number,pvalue,pathXY,saida,items+"-"+keys)        

secondIterator(8,0.05,myPathGEQC,metrics,listaGeneral)
secondIterator(8,0.05,myPathGEGH,metrics,listaGeneral)
secondIterator(11,0.05,myPathDEQC,metrics,listaDeep)
secondIterator(10,0.05,myPathDEGH,metrics,listaDeep2)

print("Finished Bonff")
#d, res = cliffs_delta(x1, x2)
def thirdIterator(path, lista, metric):

    posicao = 8

    for elements in metric:

        if elements == "lcom*":
            opt = elements.strip("*")
        else:
            opt = elements

        for items in lista:

            for keys in lista:
                if items!=keys:
                    file1 = open(path+items+".txt","r")
                    file2 = open(path+keys+".txt","r")

                    head = file1.readline()
                    head = file2.readline()
                    meuX = []
                    meuY = []
                    for lines in file1:
                        if ";NaN;" not in lines:       
                            string = lines.split(";")
                            meuX.append(float(string[posicao])+0.000000000000001)

                    for lines in file2:
                        if ";NaN;" not in lines:       
                            string = lines.split(";")
                            meuY.append(float(string[posicao])+0.000000000000001)

                    #print(path+elements+"/"+items+"-"+keys+"Cliff.txt")
                    minhaSaida = open(path+opt+"/"+items+"-"+keys+"Cliff.txt","w")

                    d, res = cliffs_delta(meuX, meuY)
                    #print(d,res)

                    minhaSaida.write(items+"-"+keys+";"+str(d)+";"+str(res))
                    
                    file1.close()
                    file2.close()
                    minhaSaida.close()

        posicao = posicao + 1

thirdIterator(myPathGEQC,listaGeneral,metrics)
thirdIterator(myPathDEQC,listaDeep,metrics)
thirdIterator(myPathGEGH,listaGeneral,metrics)
thirdIterator(myPathDEGH,listaDeep2,metrics)

print("finish cliff")

def unifyData(path,metrics,lista,tipo):

    for elements in metrics:
        if elements == "lcom*":
            opt = elements.strip("*")
        else:
            opt = elements
        saidaU = open(path+opt+"/"+tipo+"Unified.txt","w")
        for keys in lista:
            for items in lista:
                if keys!=items:
                    file = open(path+opt+"/"+keys+"-"+items+tipo+".txt","r")
                    for lines in file:
                        saidaU.write(lines+"\n")
                    file.close()

        saidaU.close()

unifyData(myPathGEQC,metrics,listaGeneral,"Cliff")
unifyData(myPathDEQC,metrics,listaDeep,"Cliff")
unifyData(myPathGEGH,metrics,listaGeneral,"Cliff")
unifyData(myPathDEGH,metrics,listaDeep2,"Cliff")


def unifyFiles(path,metric):

    arq1 = "CliffUnified.txt"
    arq2 = "Bonferroni.txt"

    for elements in metric:
        if elements == "lcom*":
            opt = elements.strip("*")
        else:
            opt = elements
        input1 = open(path+opt+"/"+arq1,"r")
        input2 = open(path+opt+"/"+arq2,"r")

        saida = open(path+opt+"/"+"BonfAndCliff.txt","w")

        linesFrom1 = []
        linesFrom2 = []
        
        for lines in input1:
            linesFrom1.append(lines.strip("\n"))

        for lines in input2:
            linesFrom2.append(lines.strip("\n"))
            

        for i in range(0,len(linesFrom1)):
            saida.write(linesFrom1[i]+";"+linesFrom2[i]+"\n")
    
        saida.close()
        input1.close()
        input2.close()
        
unifyFiles(myPathGEQC,metrics)
unifyFiles(myPathDEQC,metrics)
unifyFiles(myPathGEGH,metrics)
unifyFiles(myPathDEGH,metrics)

def turnAbsolute(path,metric):
    for elements in metric:
        if elements == "lcom*":
            opt = elements.strip("*")
        else:
            opt = elements
        entrada = open(path+opt+"/BonfAndCliff.txt","r")
        saida = open(path+opt+"/BonfAndCliffAbsolute.txt","w")        

        for lines in entrada:

            myString = lines.split(";")
            #print(myString)

            if float(myString[1])<0:
                myAbsoluteValue = float(myString[1])*(-1)
            else:
                myAbsoluteValue = float(myString[1])

            saida.write(myString[0]+";"+str(myAbsoluteValue)+";"+myString[1]+";"+myString[2]+";"
                        +myString[3]+";"+myString[4]+";"+myString[5]+";"+myString[6])            
        
        entrada.close()
        saida.close()

turnAbsolute(myPathGEQC,metrics)
turnAbsolute(myPathDEQC,metrics)
turnAbsolute(myPathGEGH,metrics)
turnAbsolute(myPathDEGH,metrics)

def addSameAggConstant(path,metric,lista):

    for elements in metric:
        if elements == "lcom*":
            opt = elements.strip("*")
        else:
            opt = elements
        entrada = open(path+opt+"/BonfAndCliffAbsolute.txt","r")
        saida = open(path+opt+"/BonfAndCliffComplete.txt","w")
        saida.write("aggType;cdAbsolute;cdOriginal;cdClassification;aggType2;mannUD;threshold\n")

        for items in lista:
            saida.write(items+";"+items+";-100;-100;invalid;"+items+"-"+items+";invalid;-100;0.005\n")
        for lines in entrada:
            saida.write(lines)
        
        entrada.close()
        saida.close()
    
addSameAggConstant(myPathGEQC,metrics,listaGeneral)
addSameAggConstant(myPathDEQC,metrics,listaDeep)
addSameAggConstant(myPathGEGH,metrics,listaGeneral)
addSameAggConstant(myPathDEGH,metrics,listaDeep2)

def tiraBarra(path,metric):

    for elements in metric:
        if elements == "lcom*":
            opt = elements.strip("*")
        else:
            opt = elements
        entrada = open(path+opt+"/BonfAndCliffComplete.txt","r")
        header = entrada.readline()
        saida = open(path+opt+"/BonfAndCliffCo.txt","w")
        saida.write("aggType,aggPar,cdAbsolute,cdOriginal,cdClassification,aggType2,manClass,mannUD,threshold\n")

        for lines in entrada:
            myList = lines.split(";")
            myString = myList[0].replace("-",",")
            final = myList[len(myList)-1].strip("\n")
            if len(myList) == 9:
                myLine = lines.replace(";",",")
                saida.write(myLine)

            else:
                saida.write(myString+","+myList[1]+","+myList[2]+","+myList[3]+","+myList[4]+","+myList[5]+","+myList[6]+","+final+"\n")

        entrada.close()
        saida.close()


tiraBarra(myPathGEQC,metrics)
tiraBarra(myPathDEQC,metrics)
tiraBarra(myPathGEGH,metrics)
tiraBarra(myPathDEGH,metrics)

def geraCSVs(path,metrics):

    for elements in metrics:
        if elements == "lcom*":
            opt = elements.strip("*")
        else:
            opt = elements
        entrada = open(path+opt+"/BonfAndCliffCo.txt","r")
        header = entrada.readline()
        saida = open(path+"myHeatMap/"+opt+"X.csv","w")
        saida.write(header)

        for lines in entrada:
            saida.write(lines)

        entrada.close()
        saida.close()

geraCSVs(myPathGEQC,metrics)
geraCSVs(myPathDEQC,metrics)
geraCSVs(myPathGEGH,metrics)
geraCSVs(myPathDEGH,metrics)

print("finishing generating csvs")

def devolveNome(nome):

    if nome=="heterogeneous":
        return "Heterogeneous"
    elif nome=="aggFE":
        return "HomogeneousFE"
    elif nome=="aggLM":
        return "HomogeneousLM"
    elif nome=="uniqueLC":
        return "IsolatedLC"
    elif nome=="uniqueRB":
        return "IsolatedRB"
    elif nome=="uniqueFE":
        return "IsolatedFE"
    elif nome=="uniqueLM":
        return "IsolatedLM"
    elif nome=="wbs":
        return "Clean"

    

def colocaNomesCertos(path,metrics):

    for elements in metrics:
        if elements == "lcom*":
            opt = elements.strip("*")
        else:
            opt = elements
        entrada = open(path+"myHeatMap/"+opt+"X.csv","r")
        saida = open(path+"myHeatMap/"+opt+".csv","w")
        header = entrada.readline()
        saida.write(header)
        for items in entrada:
            myList = items.split(",")
            #print(myList)
            name1 = devolveNome(myList[0])
            name2 = devolveNome(myList[1])
            saida.write(name1+','+name2+','+myList[2]+","+myList[3]+","+myList[4]+','+myList[5]+','+myList[6]+','+myList[7]+','+myList[8])

        entrada.close()
        saida.close()

colocaNomesCertos(myPathGEQC,metrics)
colocaNomesCertos(myPathGEGH,metrics)

print("finish details")
