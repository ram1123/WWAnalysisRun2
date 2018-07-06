def closest(list, Number):
    aux = []
    for valor in list:
        aux.append(abs(Number- float(valor)))

    return aux.index(min(aux))

#file1 = open("test.txt")
#file1 = open("moreStat.txt")
#file1 = open("Neutrino_pt.txt")
file1 = open("Neutrino_pt_type02.txt")

Total=0
CountType0 = 0
CountType1 = 0
CountType2 = 0
CountComplex = 0
#print "Truth \t Sol1 \t Sol2 \t type0 \t type1 \t type2"
#print "Truth \t Sol1 \t Sol2 \t type0 \t type2"
print "IsComplex, Truth , Sol1 , Sol2 , type0 , type1, type2, Expected soultion"
for line in file1:
   if Total == 0:	# Skip first line in txt file as it contains header
   	Total += 1
	continue;
   newLine = line.split()
   newList = newLine[2:4]	# This list contains solution of quadratic solution
   newList2 = newLine[4:]	# This contains different types of solution picked
   ClosestSol = closest(newList, float(newLine[1]))
   #print "==> ",newLine,"\t",newList2,"\t",newList[ClosestSol]
   #print "==> ",newLine,"\t == \t",newList[ClosestSol]
   print newLine[0],",",newLine[1],",",newLine[2],",",newLine[3],",",newLine[4],",",newLine[5],",",newLine[6],",",newList[ClosestSol]
   if newList[ClosestSol] == newList2[0]:
   	CountType0 += 1
   if newList[ClosestSol] == newList2[1]:
   	CountType1 += 1
   if newList[ClosestSol] == newList2[2]:
   	CountType2 += 1
   if newLine[0] == "1":
   	CountComplex += 1
	#print newLine
   Total += 1

print "% of wrong solution (type0)  = ",((Total - CountType0)/float(Total))*100.0
print "% of wrong solution (type1)  = ",((Total - CountType1)/float(Total))*100.0
print "% of wrong solution (type2)  = ",((Total - CountType2)/float(Total))*100.0
print "% of complex occurence  = ",((CountComplex)/float(Total))*100.0
