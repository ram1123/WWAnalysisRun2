def closest(list, Number):
    aux = []
    for valor in list:
        aux.append(abs(Number- float(valor)))

    return aux.index(min(aux))

#file1 = open("test.txt")
#file1 = open("moreStat.txt")
file1 = open("Neutrino_pt.txt")

Total=0
CountType0 = 0
CountType2 = 0
print "GeneratorVal   Type0 \t Type2 \t OtherSol"
for line in file1:
   if Total == 0: 
	Total += 1
   	continue
   newLine = line.split()
   newList = newLine[1:3]
   newList2 = [newLine[1]]+newLine[3:]
   ClosestSol = closest(newList, float(newLine[0]))
   print "==> ",newLine,"\t",newList,"\t",newList[ClosestSol]
   if newList[ClosestSol] == newList2[0]:
   	CountType0 += 1
   if newList[ClosestSol] == newList2[1]:
   	CountType2 += 1
   Total += 1

print "% of wrong solution (type0)  = ",((Total - CountType0)/float(Total))*100.0
print "% of wrong solution (type2)  = ",((Total - CountType2)/float(Total))*100.0
