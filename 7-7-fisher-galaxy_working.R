oldw <- getOption("warn")
options(warn = -1)
#the Galaxy tool will error out if the system displays a warning, even if that warning does not cause the code to fail
#in order to ensure that the tool does not error out due to a simple Warning, I turn down the warning option to ignore warnings

PositiveSubstrateList<- read.csv("substrates.csv", stringsAsFactors=FALSE)
NegativeSubstrateList<- read.csv("negatives.csv", stringsAsFactors=FALSE)
SubstrateBackgroundFrequency<- read.csv("SBF.csv", stringsAsFactors=FALSE,header = FALSE)
SubstrateBackgroundFrequency<-t(SubstrateBackgroundFrequency)
SubstrateBackgroundFrequency<-SubstrateBackgroundFrequency[2:nrow(SubstrateBackgroundFrequency),]
#this takes in all the inputs and ensures they are oriented correctly

ScreenerFilename<-"screener"
screaner<-read.csv(ScreenerFilename, header = FALSE, stringsAsFactors = FALSE)
#the screener file is packaged with the tool and is used for actually not for anything

DataFilename<-"thedata.RData"
load(DataFilename)
#this is a data file containing a list called Genelist
#genelist is a list of every single Uniprot-Swissprot protein in the human genome
#and those proteins are then *in-silico* digested by trypsin, and any Y-containing tryptic peptides go in the list called GeneList
#so GeneList is a list of all Y-containing tryptic peptides from the human proteome

SDtableAndPercentTable<-"output1.csv"
NormalizationScore_CharacterizationTable<-"output2.csv"
SequenceScoringAndScreening<-"output3.csv"
#these are the names of the output files that will be generated

substrates<-matrix(rep("A",times=((nrow(PositiveSubstrateList)-1)*15)),ncol = 15)
#an "empty" matrix will gets filled below.  This matrix will contain all the letters of the input Positive Substrates

for (i in 2:nrow(PositiveSubstrateList))
{
  substratemotif<-PositiveSubstrateList[i,4:18]
  substratemotif[8]<-"Y"
  #substratemotif<-paste(substratemotif,sep = "",collapse = "")
  j=i-1
  substratemotif<-unlist(substratemotif)
  substrates[j,1:15]<-substratemotif
}

substrates2<-substrates
substrates2[substrates2==""]<-"O"
#any blank spots on the substrates list are places where there was no Amino Acid, so I give them the letter O.


#I create the percent table here, which finds what percertage of each amino acid was present at each position in the substrates, from -7 to +7
if (1==1){
Column1<-substrates[,1]
Column2<-substrates[,2]
Column3<-substrates[,3]
Column4<-substrates[,4]
Column5<-substrates[,5]
Column6<-substrates[,6]
Column7<-substrates[,7]
Column8<-substrates[,8]
Column9<-substrates[,9]
Column10<-substrates[,10]
Column11<-substrates[,11]
Column12<-substrates[,12]
Column13<-substrates[,13]
Column14<-substrates[,14]
Column15<-substrates[,15]

spaces1<-sum(Column1%in% "")
spaces2<-sum(Column2%in% "")
spaces3<-sum(Column3%in% "")
spaces4<-sum(Column4%in% "")
spaces5<-sum(Column5%in% "")
spaces6<-sum(Column6%in% "")
spaces7<-sum(Column7%in% "")
spaces8<-sum(Column8%in% "")
spaces9<-sum(Column9%in% "")
spaces10<-sum(Column10%in% "")
spaces11<-sum(Column11%in% "")
spaces12<-sum(Column12%in% "")
spaces13<-sum(Column13%in% "")
spaces14<-sum(Column14%in% "")
spaces15<-sum(Column15%in% "")
OllOs<-cbind(spaces1,spaces2,spaces3,spaces4,spaces5,spaces6,spaces7,spaces8,spaces9,spaces10,spaces11,
             spaces12,spaces13,spaces14,spaces15)

A1<-sum(Column1 %in% "A")
A2<-sum(Column2 %in% "A")
A3<-sum(Column3 %in% "A")
A4<-sum(Column4 %in% "A")
A5<-sum(Column5 %in% "A")
A6<-sum(Column6 %in% "A")
A7<-sum(Column7 %in% "A")
A8<-sum(Column8 %in% "A")
A9<-sum(Column9 %in% "A")
A10<-sum(Column10 %in% "A")
A11<-sum(Column11 %in% "A")
A12<-sum(Column12 %in% "A")
A13<-sum(Column13 %in% "A")
A14<-sum(Column14 %in% "A")
A15<-sum(Column15 %in% "A")
AllAs<-cbind(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15)

C1<-sum(Column1 %in% "C")
C2<-sum(Column2 %in% "C")
C3<-sum(Column3 %in% "C")
C4<-sum(Column4 %in% "C")
C5<-sum(Column5 %in% "C")
C6<-sum(Column6 %in% "C")
C7<-sum(Column7 %in% "C")
C8<-sum(Column8 %in% "C")
C9<-sum(Column9 %in% "C")
C10<-sum(Column10 %in% "C")
C11<-sum(Column11 %in% "C")
C12<-sum(Column12 %in% "C")
C13<-sum(Column13 %in% "C")
C14<-sum(Column14 %in% "C")
C15<-sum(Column15 %in% "C")
CllCs<-cbind(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15)

D1<-sum(Column1 %in% "D")
D2<-sum(Column2 %in% "D")
D3<-sum(Column3 %in% "D")
D4<-sum(Column4 %in% "D")
D5<-sum(Column5 %in% "D")
D6<-sum(Column6 %in% "D")
D7<-sum(Column7 %in% "D")
D8<-sum(Column8 %in% "D")
D9<-sum(Column9 %in% "D")
D10<-sum(Column10 %in% "D")
D11<-sum(Column11 %in% "D")
D12<-sum(Column12 %in% "D")
D13<-sum(Column13 %in% "D")
D14<-sum(Column14 %in% "D")
D15<-sum(Column15 %in% "D")
DllDs<-cbind(D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15)

E1<-sum(Column1 %in% "E")
E2<-sum(Column2 %in% "E")
E3<-sum(Column3 %in% "E")
E4<-sum(Column4 %in% "E")
E5<-sum(Column5 %in% "E")
E6<-sum(Column6 %in% "E")
E7<-sum(Column7 %in% "E")
E8<-sum(Column8 %in% "E")
E9<-sum(Column9 %in% "E")
E10<-sum(Column10 %in% "E")
E11<-sum(Column11 %in% "E")
E12<-sum(Column12 %in% "E")
E13<-sum(Column13 %in% "E")
E14<-sum(Column14 %in% "E")
E15<-sum(Column15 %in% "E")
EllEs<-cbind(E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13,E14,E15)

F1<-sum(Column1 %in% "F")
F2<-sum(Column2 %in% "F")
F3<-sum(Column3 %in% "F")
F4<-sum(Column4 %in% "F")
F5<-sum(Column5 %in% "F")
F6<-sum(Column6 %in% "F")
F7<-sum(Column7 %in% "F")
F8<-sum(Column8 %in% "F")
F9<-sum(Column9 %in% "F")
F10<-sum(Column10 %in% "F")
F11<-sum(Column11 %in% "F")
F12<-sum(Column12 %in% "F")
F13<-sum(Column13 %in% "F")
F14<-sum(Column14 %in% "F")
F15<-sum(Column15 %in% "F")
FllFs<-cbind(F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15)

G1<-sum(Column1 %in% "G")
G2<-sum(Column2 %in% "G")
G3<-sum(Column3 %in% "G")
G4<-sum(Column4 %in% "G")
G5<-sum(Column5 %in% "G")
G6<-sum(Column6 %in% "G")
G7<-sum(Column7 %in% "G")
G8<-sum(Column8 %in% "G")
G9<-sum(Column9 %in% "G")
G10<-sum(Column10 %in% "G")
G11<-sum(Column11 %in% "G")
G12<-sum(Column12 %in% "G")
G13<-sum(Column13 %in% "G")
G14<-sum(Column14 %in% "G")
G15<-sum(Column15 %in% "G")
GllGs<-cbind(G1,G2,G3,G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G14,G15)

H1<-sum(Column1 %in% "H")
H2<-sum(Column2 %in% "H")
H3<-sum(Column3 %in% "H")
H4<-sum(Column4 %in% "H")
H5<-sum(Column5 %in% "H")
H6<-sum(Column6 %in% "H")
H7<-sum(Column7 %in% "H")
H8<-sum(Column8 %in% "H")
H9<-sum(Column9 %in% "H")
H10<-sum(Column10 %in% "H")
H11<-sum(Column11 %in% "H")
H12<-sum(Column12 %in% "H")
H13<-sum(Column13 %in% "H")
H14<-sum(Column14 %in% "H")
H15<-sum(Column15 %in% "H")
HllHs<-cbind(H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,H11,H12,H13,H14,H15)

I1<-sum(Column1 %in% "I")
I2<-sum(Column2 %in% "I")
I3<-sum(Column3 %in% "I")
I4<-sum(Column4 %in% "I")
I5<-sum(Column5 %in% "I")
I6<-sum(Column6 %in% "I")
I7<-sum(Column7 %in% "I")
I8<-sum(Column8 %in% "I")
I9<-sum(Column9 %in% "I")
I10<-sum(Column10 %in% "I")
I11<-sum(Column11 %in% "I")
I12<-sum(Column12 %in% "I")
I13<-sum(Column13 %in% "I")
I14<-sum(Column14 %in% "I")
I15<-sum(Column15 %in% "I")
IllIs<-cbind(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I13,I14,I15)

K1<-sum(Column1 %in% "K")
K2<-sum(Column2 %in% "K")
K3<-sum(Column3 %in% "K")
K4<-sum(Column4 %in% "K")
K5<-sum(Column5 %in% "K")
K6<-sum(Column6 %in% "K")
K7<-sum(Column7 %in% "K")
K8<-sum(Column8 %in% "K")
K9<-sum(Column9 %in% "K")
K10<-sum(Column10 %in% "K")
K11<-sum(Column11 %in% "K")
K12<-sum(Column12 %in% "K")
K13<-sum(Column13 %in% "K")
K14<-sum(Column14 %in% "K")
K15<-sum(Column15 %in% "K")
KllKs<-cbind(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15)

L1<-sum(Column1 %in% "L")
L2<-sum(Column2 %in% "L")
L3<-sum(Column3 %in% "L")
L4<-sum(Column4 %in% "L")
L5<-sum(Column5 %in% "L")
L6<-sum(Column6 %in% "L")
L7<-sum(Column7 %in% "L")
L8<-sum(Column8 %in% "L")
L9<-sum(Column9 %in% "L")
L10<-sum(Column10 %in% "L")
L11<-sum(Column11 %in% "L")
L12<-sum(Column12 %in% "L")
L13<-sum(Column13 %in% "L")
L14<-sum(Column14 %in% "L")
L15<-sum(Column15 %in% "L")
LllLs<-cbind(L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15)

M1<-sum(Column1 %in% "M")
M2<-sum(Column2 %in% "M")
M3<-sum(Column3 %in% "M")
M4<-sum(Column4 %in% "M")
M5<-sum(Column5 %in% "M")
M6<-sum(Column6 %in% "M")
M7<-sum(Column7 %in% "M")
M8<-sum(Column8 %in% "M")
M9<-sum(Column9 %in% "M")
M10<-sum(Column10 %in% "M")
M11<-sum(Column11 %in% "M")
M12<-sum(Column12 %in% "M")
M13<-sum(Column13 %in% "M")
M14<-sum(Column14 %in% "M")
M15<-sum(Column15 %in% "M")
MllMs<-cbind(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,M15)

N1<-sum(Column1 %in% "N")
N2<-sum(Column2 %in% "N")
N3<-sum(Column3 %in% "N")
N4<-sum(Column4 %in% "N")
N5<-sum(Column5 %in% "N")
N6<-sum(Column6 %in% "N")
N7<-sum(Column7 %in% "N")
N8<-sum(Column8 %in% "N")
N9<-sum(Column9 %in% "N")
N10<-sum(Column10 %in% "N")
N11<-sum(Column11 %in% "N")
N12<-sum(Column12 %in% "N")
N13<-sum(Column13 %in% "N")
N14<-sum(Column14 %in% "N")
N15<-sum(Column15 %in% "N")
NllNs<-cbind(N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15)

P1<-sum(Column1 %in% "P")
P2<-sum(Column2 %in% "P")
P3<-sum(Column3 %in% "P")
P4<-sum(Column4 %in% "P")
P5<-sum(Column5 %in% "P")
P6<-sum(Column6 %in% "P")
P7<-sum(Column7 %in% "P")
P8<-sum(Column8 %in% "P")
P9<-sum(Column9 %in% "P")
P10<-sum(Column10 %in% "P")
P11<-sum(Column11 %in% "P")
P12<-sum(Column12 %in% "P")
P13<-sum(Column13 %in% "P")
P14<-sum(Column14 %in% "P")
P15<-sum(Column15 %in% "P")
PllPs<-cbind(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15)

Q1<-sum(Column1 %in% "Q")
Q2<-sum(Column2 %in% "Q")
Q3<-sum(Column3 %in% "Q")
Q4<-sum(Column4 %in% "Q")
Q5<-sum(Column5 %in% "Q")
Q6<-sum(Column6 %in% "Q")
Q7<-sum(Column7 %in% "Q")
Q8<-sum(Column8 %in% "Q")
Q9<-sum(Column9 %in% "Q")
Q10<-sum(Column10 %in% "Q")
Q11<-sum(Column11 %in% "Q")
Q12<-sum(Column12 %in% "Q")
Q13<-sum(Column13 %in% "Q")
Q14<-sum(Column14 %in% "Q")
Q15<-sum(Column15 %in% "Q")
QllQs<-cbind(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15)

R1<-sum(Column1 %in% "R")
R2<-sum(Column2 %in% "R")
R3<-sum(Column3 %in% "R")
R4<-sum(Column4 %in% "R")
R5<-sum(Column5 %in% "R")
R6<-sum(Column6 %in% "R")
R7<-sum(Column7 %in% "R")
R8<-sum(Column8 %in% "R")
R9<-sum(Column9 %in% "R")
R10<-sum(Column10 %in% "R")
R11<-sum(Column11 %in% "R")
R12<-sum(Column12 %in% "R")
R13<-sum(Column13 %in% "R")
R14<-sum(Column14 %in% "R")
R15<-sum(Column15 %in% "R")
RllRs<-cbind(R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15)

S1<-sum(Column1 %in% "S")
S2<-sum(Column2 %in% "S")
S3<-sum(Column3 %in% "S")
S4<-sum(Column4 %in% "S")
S5<-sum(Column5 %in% "S")
S6<-sum(Column6 %in% "S")
S7<-sum(Column7 %in% "S")
S8<-sum(Column8 %in% "S")
S9<-sum(Column9 %in% "S")
S10<-sum(Column10 %in% "S")
S11<-sum(Column11 %in% "S")
S12<-sum(Column12 %in% "S")
S13<-sum(Column13 %in% "S")
S14<-sum(Column14 %in% "S")
S15<-sum(Column15 %in% "S")
SllSs<-cbind(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15)

T1<-sum(Column1 %in% "T")
T2<-sum(Column2 %in% "T")
T3<-sum(Column3 %in% "T")
T4<-sum(Column4 %in% "T")
T5<-sum(Column5 %in% "T")
T6<-sum(Column6 %in% "T")
T7<-sum(Column7 %in% "T")
T8<-sum(Column8 %in% "T")
T9<-sum(Column9 %in% "T")
T10<-sum(Column10 %in% "T")
T11<-sum(Column11 %in% "T")
T12<-sum(Column12 %in% "T")
T13<-sum(Column13 %in% "T")
T14<-sum(Column14 %in% "T")
T15<-sum(Column15 %in% "T")
TllTs<-cbind(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15)

V1<-sum(Column1 %in% "V")
V2<-sum(Column2 %in% "V")
V3<-sum(Column3 %in% "V")
V4<-sum(Column4 %in% "V")
V5<-sum(Column5 %in% "V")
V6<-sum(Column6 %in% "V")
V7<-sum(Column7 %in% "V")
V8<-sum(Column8 %in% "V")
V9<-sum(Column9 %in% "V")
V10<-sum(Column10 %in% "V")
V11<-sum(Column11 %in% "V")
V12<-sum(Column12 %in% "V")
V13<-sum(Column13 %in% "V")
V14<-sum(Column14 %in% "V")
V15<-sum(Column15 %in% "V")
VllVs<-cbind(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15)

W1<-sum(Column1 %in% "W")
W2<-sum(Column2 %in% "W")
W3<-sum(Column3 %in% "W")
W4<-sum(Column4 %in% "W")
W5<-sum(Column5 %in% "W")
W6<-sum(Column6 %in% "W")
W7<-sum(Column7 %in% "W")
W8<-sum(Column8 %in% "W")
W9<-sum(Column9 %in% "W")
W10<-sum(Column10 %in% "W")
W11<-sum(Column11 %in% "W")
W12<-sum(Column12 %in% "W")
W13<-sum(Column13 %in% "W")
W14<-sum(Column14 %in% "W")
W15<-sum(Column15 %in% "W")
WllWs<-cbind(W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,W13,W14,W15)

Y1<-sum(Column1 %in% "Y")
Y2<-sum(Column2 %in% "Y")
Y3<-sum(Column3 %in% "Y")
Y4<-sum(Column4 %in% "Y")
Y5<-sum(Column5 %in% "Y")
Y6<-sum(Column6 %in% "Y")
Y7<-sum(Column7 %in% "Y")
Y8<-sum(Column8 %in% "Y")
Y9<-sum(Column9 %in% "Y")
Y10<-sum(Column10 %in% "Y")
Y11<-sum(Column11 %in% "Y")
Y12<-sum(Column12 %in% "Y")
Y13<-sum(Column13 %in% "Y")
Y14<-sum(Column14 %in% "Y")
Y15<-sum(Column15 %in% "Y")
YllYs<-cbind(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,Y11,Y12,Y13,Y14,Y15)

PercentTable<-rbind(AllAs,CllCs,DllDs,EllEs,FllFs,GllGs,HllHs,IllIs,KllKs,LllLs,MllMs,NllNs,PllPs,QllQs,RllRs,SllSs,TllTs,VllVs,WllWs,YllYs,OllOs)
}

AllSubBackFreq<-array(data = NA,dim = c(21,15,nrow(SubstrateBackgroundFrequency)))
#this is where I'm creating the new SubBackFreq table, I have a list with all possible SBF matrices,
#I perform a for function to find all the matrices that are in Substrate Background Frequency
#and place them all in this array, then I will do mean and SD

AAccessionNumbers<-SubstrateBackgroundFrequency[,1]
AllGeneNames<-names(Genelist)

number_replaced<-0
totalmotifs<-0
for (z in 1:length(AAccessionNumbers)) {
  #for every protein found in the substrate background frequency list
  pattern<-AAccessionNumbers[z]
  referencepoint<-grepl(pattern = pattern, x=AllGeneNames,fixed = TRUE)
  #find that protein in the GeneList list
  referencenumber<-which(referencepoint==TRUE)
  if (length(referencenumber)<1){referencenumber<-FALSE}
  if (referencenumber!=FALSE){
    #then if you have actually found a match in the GeneList list, you will do the following
    motifs<-unlist(Genelist[[referencenumber]])
    #take every tryptic peptide from that list
    therow<-c(1:15)
    for (a in 1:length(motifs)) {
      thecut<-unlist(strsplit(motifs[a], split=""))
      edges<-c("O","O","O","O","O","O","O")
      thecut<-c(edges,thecut,edges)
      theYs<-which(thecut=="Y")
      #and use the above code to replace blank spaces in the -7 to +7 motif with Os, same as what I did for the positive substrates
      for (q in 1:length(theYs)) {
        thiscut<-thecut[(theYs[q]-7):(theYs[q]+7)]
        therow<-rbind(therow,thiscut)
        totalmotifs<-totalmotifs+1
        #count how many motifs you've put in there
      }
    }
    
    cutreplacement<-c("X","X","X","X","X","X","X","X","X","X","X","X","X","X","X")
    for (t in 1:nrow(therow)) {
      compare1<-therow[t,1:15]
      compare1<-paste(compare1,sep = "",collapse = "")
      
      for (v in 1:nrow(substrates2)) {
        positivesubstrate<-substrates2[v,1:15]
        positivesubstrate<-paste(positivesubstrate,sep = "",collapse = "")
        
        if (compare1==positivesubstrate){
          therow[t,1:15]<-cutreplacement
          number_replaced<-number_replaced+1
          #if any of the substrates found in the substrate background frequency are identical to the positive substrates, this removes them
          #and replaces them with Xs to ensure they do not get counted when all the AAs are being counted
          }
      }
      
    }
    
    Column1<-therow[,1]
    Column2<-therow[,2]
    Column3<-therow[,3]
    Column4<-therow[,4]
    Column5<-therow[,5]
    Column6<-therow[,6]
    Column7<-therow[,7]
    Column8<-therow[,8]
    Column9<-therow[,9]
    Column10<-therow[,10]
    Column11<-therow[,11]
    Column12<-therow[,12]
    Column13<-therow[,13]
    Column14<-therow[,14]
    Column15<-therow[,15]
    #split the motifs you found above into individual columns
    #each column is ome of the positions between -7 and +7, and so each column represents every single amino acid found at that position
    slice1<-c(sum(Column1=="A"),sum(Column1=="C"),sum(Column1=="D"),sum(Column1=="E"),sum(Column1=="F"),
              sum(Column1=="G"),sum(Column1=="H"),sum(Column1=="I"),sum(Column1=="K"),sum(Column1=="L"),
              sum(Column1=="M"),sum(Column1=="N"),sum(Column1=="P"),sum(Column1=="Q"),sum(Column1=="R"),
              sum(Column1=="S"),sum(Column1=="T"),sum(Column1=="V"),sum(Column1=="W"),sum(Column1=="Y"),
              sum(Column1=="O"))
    #count how many of each amino acid were present in that column, the vector of "how many of each amino acid were present in Column1", the name of that vector is slice1
    slice2<-c(sum(Column2=="A"),sum(Column2=="C"),sum(Column2=="D"),sum(Column2=="E"),sum(Column2=="F"),
              sum(Column2=="G"),sum(Column2=="H"),sum(Column2=="I"),sum(Column2=="K"),sum(Column2=="L"),
              sum(Column2=="M"),sum(Column2=="N"),sum(Column2=="P"),sum(Column2=="Q"),sum(Column2=="R"),
              sum(Column2=="S"),sum(Column2=="T"),sum(Column2=="V"),sum(Column2=="W"),sum(Column2=="Y"),
              sum(Column2=="O"))
    slice3<-c(sum(Column3=="A"),sum(Column3=="C"),sum(Column3=="D"),sum(Column3=="E"),sum(Column3=="F"),
              sum(Column3=="G"),sum(Column3=="H"),sum(Column3=="I"),sum(Column3=="K"),sum(Column3=="L"),
              sum(Column3=="M"),sum(Column3=="N"),sum(Column3=="P"),sum(Column3=="Q"),sum(Column3=="R"),
              sum(Column3=="S"),sum(Column3=="T"),sum(Column3=="V"),sum(Column3=="W"),sum(Column3=="Y"),
              sum(Column3=="O"))
    slice4<-c(sum(Column4=="A"),sum(Column4=="C"),sum(Column4=="D"),sum(Column4=="E"),sum(Column4=="F"),
              sum(Column4=="G"),sum(Column4=="H"),sum(Column4=="I"),sum(Column4=="K"),sum(Column4=="L"),
              sum(Column4=="M"),sum(Column4=="N"),sum(Column4=="P"),sum(Column4=="Q"),sum(Column4=="R"),
              sum(Column4=="S"),sum(Column4=="T"),sum(Column4=="V"),sum(Column4=="W"),sum(Column4=="Y"),
              sum(Column4=="O"))
    slice5<-c(sum(Column5=="A"),sum(Column5=="C"),sum(Column5=="D"),sum(Column5=="E"),sum(Column5=="F"),
              sum(Column5=="G"),sum(Column5=="H"),sum(Column5=="I"),sum(Column5=="K"),sum(Column5=="L"),
              sum(Column5=="M"),sum(Column5=="N"),sum(Column5=="P"),sum(Column5=="Q"),sum(Column5=="R"),
              sum(Column5=="S"),sum(Column5=="T"),sum(Column5=="V"),sum(Column5=="W"),sum(Column5=="Y"),
              sum(Column5=="O"))
    slice6<-c(sum(Column6=="A"),sum(Column6=="C"),sum(Column6=="D"),sum(Column6=="E"),sum(Column6=="F"),
              sum(Column6=="G"),sum(Column6=="H"),sum(Column6=="I"),sum(Column6=="K"),sum(Column6=="L"),
              sum(Column6=="M"),sum(Column6=="N"),sum(Column6=="P"),sum(Column6=="Q"),sum(Column6=="R"),
              sum(Column6=="S"),sum(Column6=="T"),sum(Column6=="V"),sum(Column6=="W"),sum(Column6=="Y"),
              sum(Column6=="O"))
    slice7<-c(sum(Column7=="A"),sum(Column7=="C"),sum(Column7=="D"),sum(Column7=="E"),sum(Column7=="F"),
              sum(Column7=="G"),sum(Column7=="H"),sum(Column7=="I"),sum(Column7=="K"),sum(Column7=="L"),
              sum(Column7=="M"),sum(Column7=="N"),sum(Column7=="P"),sum(Column7=="Q"),sum(Column7=="R"),
              sum(Column7=="S"),sum(Column7=="T"),sum(Column7=="V"),sum(Column7=="W"),sum(Column7=="Y"),
              sum(Column7=="O"))
    slice8<-c(sum(Column8=="A"),sum(Column8=="C"),sum(Column8=="D"),sum(Column8=="E"),sum(Column8=="F"),
              sum(Column8=="G"),sum(Column8=="H"),sum(Column8=="I"),sum(Column8=="K"),sum(Column8=="L"),
              sum(Column8=="M"),sum(Column8=="N"),sum(Column8=="P"),sum(Column8=="Q"),sum(Column8=="R"),
              sum(Column8=="S"),sum(Column8=="T"),sum(Column8=="V"),sum(Column8=="W"),sum(Column8=="Y"),
              sum(Column8=="O"))
    slice9<-c(sum(Column9=="A"),sum(Column9=="C"),sum(Column9=="D"),sum(Column9=="E"),sum(Column9=="F"),
              sum(Column9=="G"),sum(Column9=="H"),sum(Column9=="I"),sum(Column9=="K"),sum(Column9=="L"),
              sum(Column9=="M"),sum(Column9=="N"),sum(Column9=="P"),sum(Column9=="Q"),sum(Column9=="R"),
              sum(Column9=="S"),sum(Column9=="T"),sum(Column9=="V"),sum(Column9=="W"),sum(Column9=="Y"),
              sum(Column9=="O"))
    slice10<-c(sum(Column10=="A"),sum(Column10=="C"),sum(Column10=="D"),sum(Column10=="E"),sum(Column10=="F"),
              sum(Column10=="G"),sum(Column10=="H"),sum(Column10=="I"),sum(Column10=="K"),sum(Column10=="L"),
              sum(Column10=="M"),sum(Column10=="N"),sum(Column10=="P"),sum(Column10=="Q"),sum(Column10=="R"),
              sum(Column10=="S"),sum(Column10=="T"),sum(Column10=="V"),sum(Column10=="W"),sum(Column10=="Y"),
              sum(Column10=="O"))
    slice11<-c(sum(Column11=="A"),sum(Column11=="C"),sum(Column11=="D"),sum(Column11=="E"),sum(Column11=="F"),
              sum(Column11=="G"),sum(Column11=="H"),sum(Column11=="I"),sum(Column11=="K"),sum(Column11=="L"),
              sum(Column11=="M"),sum(Column11=="N"),sum(Column11=="P"),sum(Column11=="Q"),sum(Column11=="R"),
              sum(Column11=="S"),sum(Column11=="T"),sum(Column11=="V"),sum(Column11=="W"),sum(Column11=="Y"),
              sum(Column11=="O"))
    slice12<-c(sum(Column12=="A"),sum(Column12=="C"),sum(Column12=="D"),sum(Column12=="E"),sum(Column12=="F"),
              sum(Column12=="G"),sum(Column12=="H"),sum(Column12=="I"),sum(Column12=="K"),sum(Column12=="L"),
              sum(Column12=="M"),sum(Column12=="N"),sum(Column12=="P"),sum(Column12=="Q"),sum(Column12=="R"),
              sum(Column12=="S"),sum(Column12=="T"),sum(Column12=="V"),sum(Column12=="W"),sum(Column12=="Y"),
              sum(Column12=="O"))
    slice13<-c(sum(Column13=="A"),sum(Column13=="C"),sum(Column13=="D"),sum(Column13=="E"),sum(Column13=="F"),
              sum(Column13=="G"),sum(Column13=="H"),sum(Column13=="I"),sum(Column13=="K"),sum(Column13=="L"),
              sum(Column13=="M"),sum(Column13=="N"),sum(Column13=="P"),sum(Column13=="Q"),sum(Column13=="R"),
              sum(Column13=="S"),sum(Column13=="T"),sum(Column13=="V"),sum(Column13=="W"),sum(Column13=="Y"),
              sum(Column13=="O"))
    slice14<-c(sum(Column14=="A"),sum(Column14=="C"),sum(Column14=="D"),sum(Column14=="E"),sum(Column14=="F"),
              sum(Column14=="G"),sum(Column14=="H"),sum(Column14=="I"),sum(Column14=="K"),sum(Column14=="L"),
              sum(Column14=="M"),sum(Column14=="N"),sum(Column14=="P"),sum(Column14=="Q"),sum(Column14=="R"),
              sum(Column14=="S"),sum(Column14=="T"),sum(Column14=="V"),sum(Column14=="W"),sum(Column14=="Y"),
              sum(Column14=="O"))
    slice15<-c(sum(Column15=="A"),sum(Column15=="C"),sum(Column15=="D"),sum(Column15=="E"),sum(Column15=="F"),
              sum(Column15=="G"),sum(Column15=="H"),sum(Column15=="I"),sum(Column15=="K"),sum(Column15=="L"),
              sum(Column15=="M"),sum(Column15=="N"),sum(Column15=="P"),sum(Column15=="Q"),sum(Column15=="R"),
              sum(Column15=="S"),sum(Column15=="T"),sum(Column15=="V"),sum(Column15=="W"),sum(Column15=="Y"),
              sum(Column15=="O"))
    ThisMatix<-cbind(slice1,slice2,slice3,slice4,slice5,slice6,slice7,slice8,slice9,
                     slice10,slice11,slice12,slice13,slice14,slice15)
    ThisMatix<-ThisMatix
    #stick all those slices together in a matrix
    AllSubBackFreq[1:21,1:15,z]<-ThisMatix
    #add that matrix to the 3d matrix of substratebackgroundfrequency
    #the full matrix is a 3d matrix.  The Z dimension represents each Accession Number aka protein, the Y dimension represents 
    #each of the 20 amino acids plus O which is AA number 21, and the X dimension represents the positions -7 to +7
  }
}

theletters<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","O")
#these are the 20 AAs of the code plus O which represents "no amino acid present here"

SumSBF<-apply(AllSubBackFreq, c(1,2), sum, na.rm=TRUE)
#the 3d matrix above is summed across its Z dimension


fisheroddstable<-matrix(data = 1,nrow = 21,ncol = 15)
fisherpvalstable<-matrix(data = 1,nrow = 21,ncol = 15)
fisherpvalstableadjusted<-matrix(data = 1,nrow = 21,ncol = 15)
#pre-allocate vectors

for (rowas in 1:21) {
  #for each of the 20 AAs plus O
  for (colams in 1:15) {
    #for each position -7 to +7
    fishermatrix<-matrix(data=c(PercentTable[rowas,colams],nrow(substrates),SumSBF[rowas,colams],(totalmotifs-number_replaced)),nrow = 2)
    thetest<-fisher.test(x=fishermatrix)
    #perform a Fisher Test to see if the frequency of that AA in the substrate background frequency is different than the frequency in the Positive Substrates
    #remember that the frequency in the Positive Substrates is determined by the percent table
    fisheroddstable[rowas,colams]<-thetest$estimate
    fisherpvalstable[rowas,colams]<-thetest$p.value
    fisherpvalstableadjusted[rowas,colams]<-p.adjust(p=thetest$p.value,method = "fdr",n=21*15)
  }
}

fisheroddstable<-cbind.data.frame(theletters,fisheroddstable)
fisherpvalstable<-cbind.data.frame(theletters,fisherpvalstable)
fisherpvalstableadjusted<-cbind.data.frame(theletters,fisherpvalstableadjusted)

fisherupdown<-fisheroddstable

for (x in 1:21) {
  for (y in 2:16) {
    theval<-1
    testval<-fisheroddstable[x,y]
    testp<-fisherpvalstable[x,y]
    if (testp<.05){
      theval<-testval
    }
    fisherupdown[x,y]<-theval
    #determine which of the Fisher Odds found were actually significant
  }
}

write.table(x="Fisher Odds, only significant ones",file = SDtableAndPercentTable, append = TRUE,sep = ",",col.names = FALSE,row.names = FALSE)
write.table(x=fisherupdown,file = SDtableAndPercentTable, append = TRUE,sep = ",",col.names = FALSE,row.names = FALSE)
write.table(x="Fisher Odds",file = SDtableAndPercentTable, append = TRUE,sep = ",",col.names = FALSE,row.names = FALSE)
write.table(x=fisheroddstable,file = SDtableAndPercentTable, append = TRUE,sep = ",",col.names = FALSE,row.names = FALSE)
write.table(x="Fisher p.values",file = SDtableAndPercentTable, append = TRUE,sep = ",",col.names = FALSE,row.names = FALSE)
write.table(x=fisherpvalstable,file = SDtableAndPercentTable, append = TRUE,sep = ",",col.names = FALSE,row.names = FALSE)
write.table(x="Fisher p.values adjusted",file = SDtableAndPercentTable, append = TRUE,sep = ",",col.names = FALSE,row.names = FALSE)
write.table(x=fisherpvalstableadjusted,file = SDtableAndPercentTable, append = TRUE,sep = ",",col.names = FALSE,row.names = FALSE)

SetOfAAs<-c("Letter","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
SetOfAAs<-matrix(data = SetOfAAs,ncol = 1)
numberofY<-as.numeric(SubstrateBackgroundFrequency[,34])
numberofY<-numberofY[!is.na(numberofY)]
numberofPY<-as.numeric(SubstrateBackgroundFrequency[,35])
numberofPY<-numberofPY[!is.na(numberofPY)]
NormalizationScore<-sum(numberofPY)/sum(numberofY)
NormalizationScore<-c("Normalization Score",NormalizationScore)
write.table(NormalizationScore, file = SiteSelectivityTable_EndogenousProbabilityMatrix_NormalizationScore_CharacterizationTable, append = TRUE,sep = ",",row.names = FALSE, col.names = FALSE)
#output files


bareSDs<-fisherupdown[1:20,2:16]
bareSDs[20,8]<-3
bareSDs[3:4,2]<-1
bareSDs[3:4,4]<-1
bareSDs[3:4,9]<-1
bareSDs[3:4,13:14]<-1
#take the Fisher Odds table but only that those odds, remove the labels that I put of them

bareSDs[20,8]<-1



A=1
C=2
D=3
E=4
F=5
G=6
H=7
I=8
K=9
L=10
M=11
N=12
P=13
Q=14
R=15
S=16
T=17
V=18
W=19
Y=20

aa_props <- c("A"=A, "C"=C, "D"=D, "E"=E, "F"=F,"G"=G,"H"=H,"I"=I,"K"=K,"L"=L,"M"=M,"N"=N,"P"=P,"Q"=Q,"R"=R,
              "S"=S,"T"=T,"V"=V,"W"=W,"Y"=Y,"xY"=Y,"O"=21)
#this function converts letters to a number based on the equalities that I have written above
#so aa_props(c(A,C,D)) will give c(1,2,3) as an answer

ThisKinTable<-fisheroddstable

NegativeScores<-rep(NA,times=nrow(NegativeSubstrateList))
NegativeWeirdScores<-rep(NA,times=nrow(NegativeSubstrateList))

#here I am going to take all of the Negative motifs and score them, according to the Fisher Odds Tables and the statistics of KINATEST-ID
for (v in 1:nrow(NegativeSubstrateList)) {
  motif<-NegativeSubstrateList[v,2]
  motif<-unlist(strsplit(motif,""))
  motif<- gsub(" ","O",motif)  
  motif <- sapply(motif, function (x) aa_props[x])
  Scoringpeptide<-motif
  Scoringpeptide<-Scoringpeptide
  ThisKinTableScore<-as.numeric(ThisKinTable[Scoringpeptide[1],2])*ThisKinTable[as.numeric(Scoringpeptide[2]),3]*ThisKinTable[as.numeric(Scoringpeptide[3]),4]*
    ThisKinTable[as.numeric(Scoringpeptide[4]),5]*ThisKinTable[as.numeric(Scoringpeptide[5]),6]*ThisKinTable[as.numeric(Scoringpeptide[6]),7]*
    ThisKinTable[as.numeric(Scoringpeptide[7]),8]*
    #ThisKinTable[as.numeric(Scoringpeptide[8]),10]*
    ThisKinTable[as.numeric(Scoringpeptide[9]),10]*ThisKinTable[as.numeric(Scoringpeptide[10]),11]*ThisKinTable[as.numeric(Scoringpeptide[11]),12]*
    ThisKinTable[as.numeric(Scoringpeptide[12]),13]*ThisKinTable[as.numeric(Scoringpeptide[13]),14]*ThisKinTable[as.numeric(Scoringpeptide[14]),15]*
    ThisKinTable[as.numeric(Scoringpeptide[15]),16]
  NegativeScores[v]<-ThisKinTableScore
  ThisKinTableScore<-(ThisKinTableScore/(ThisKinTableScore+1/as.numeric(NormalizationScore[2])))
  NegativeWeirdScores[v]<-ThisKinTableScore*100
}

negativesubstrates<-NegativeSubstrateList[,2]
NegativeWithScores<-cbind(negativesubstrates,as.character(NegativeScores),as.character(NegativeWeirdScores))

PositiveScores<-rep(NA,times=nrow(PositiveSubstrateList))
PositiveWeirdScores<-rep(NA,times=nrow(PositiveSubstrateList))

#now I'm going to score all of the positive motifs, according to the same statistics
for (v in 1:nrow(PositiveSubstrateList)) {
  motif<-PositiveSubstrateList[v,4:18]
  motif<-unlist(motif)
  motif<- gsub("^$","O",motif)
  motif <- sapply(motif, function (x) aa_props[x])
  Scoringpeptide<-motif
  Scoringpeptide<-Scoringpeptide
  ThisKinTableScore<-as.numeric(ThisKinTable[Scoringpeptide[1],2])*ThisKinTable[as.numeric(Scoringpeptide[2]),3]*ThisKinTable[as.numeric(Scoringpeptide[3]),4]*
    ThisKinTable[as.numeric(Scoringpeptide[4]),5]*ThisKinTable[as.numeric(Scoringpeptide[5]),6]*ThisKinTable[as.numeric(Scoringpeptide[6]),7]*
    ThisKinTable[as.numeric(Scoringpeptide[7]),8]*
    #ThisKinTable[as.numeric(Scoringpeptide[8]),10]*
    ThisKinTable[as.numeric(Scoringpeptide[9]),10]*ThisKinTable[as.numeric(Scoringpeptide[10]),11]*ThisKinTable[as.numeric(Scoringpeptide[11]),12]*
    ThisKinTable[as.numeric(Scoringpeptide[12]),13]*ThisKinTable[as.numeric(Scoringpeptide[13]),14]*ThisKinTable[as.numeric(Scoringpeptide[14]),15]*
    ThisKinTable[as.numeric(Scoringpeptide[15]),16]

  PositiveScores[v]<-ThisKinTableScore
  ThisKinTableScore<-(ThisKinTableScore/(ThisKinTableScore+1/as.numeric(NormalizationScore[2])))
  PositiveWeirdScores[v]<-ThisKinTableScore*100
}

positivesubstrates<-PositiveSubstrateList[,4:18]
positivewithscores<-cbind.data.frame(positivesubstrates,PositiveScores,PositiveWeirdScores)



SetOfAAs<-c("Letter","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
SumOfSigmaAAs<-c(1:15)

for (i in 1:15){
  SumOfSigmasValue<-0
  for (j in 1:20){
    value<-0
    if (bareSDs[j,i]>1){
      k<-j+1
      value<-sum(substrates[,i]==SetOfAAs[k])
    }
    SumOfSigmasValue<-SumOfSigmasValue+value
  }
  SumOfSigmaAAs[i]<-SumOfSigmasValue
}
#this creates the Sum Of Sigmas value, also from the KINATEST-ID paper

threshold<-c(1:100,(1:9)/10,(1:9)/100,0,-.1)
threshold<-threshold[order(threshold,decreasing = TRUE)]
threshold

Truepositives<-c(1:120)
Falsenegatives<-c(1:120)
Sensitivity<-c(1:120)
TrueNegatives<-c(1:120)
FalsePositives<-c(1:120)
One_Minus_Specificity<-c(1:120)
Accuracy<-c(1:120)
MCC<-c(1:120)
EER<-c(1:120)
Precision<-c(1:120)
F_One_Half<-c(1:120)
F_One<-c(1:120)
F_Two<-c(1:120)
FalsePositiveRate<-c(1:120)

#here I am creating the Matthew's Correlation Coefficient table for each posible threshold value that could be used to separate negative motifs from positive motifs
#these statistics allow us to determine the best threshold value
for (z in 1:120) {
  thres<-threshold[z]
  Truepositives[z]<-length(PositiveWeirdScores[PositiveWeirdScores>=(thres)])
  Falsenegatives[z]<-nrow(positivesubstrates)-Truepositives[z]
  Sensitivity[z]<-Truepositives[z]/(Falsenegatives[z]+Truepositives[z])
  TrueNegatives[z]<-length(NegativeWeirdScores[NegativeWeirdScores<(thres)])
  # at thresh 100 this should be 0, because it is total minus true negatives
  FalsePositives[z]<-nrow(NegativeSubstrateList)-TrueNegatives[z]
  One_Minus_Specificity[z]<-1-(TrueNegatives[z]/(FalsePositives[z]+TrueNegatives[z]))
  Accuracy[z]<-100*(Truepositives[z]+TrueNegatives[z])/(Falsenegatives[z]+FalsePositives[z]+TrueNegatives[z]+Truepositives[z])
  MCC[z]<-((Truepositives[z]*TrueNegatives[z])-(Falsenegatives[z]*FalsePositives[z]))/sqrt(round(round(Truepositives[z]+Falsenegatives[z])*round(TrueNegatives[z]+FalsePositives[z])*round(Truepositives[z]+FalsePositives[z])*round(TrueNegatives[z]+Falsenegatives[z])))
  EER[z]<-(FalsePositives[z]+Falsenegatives[z])/(Truepositives[z]+TrueNegatives[z]+FalsePositives[z]+Falsenegatives[z])
  Precision[z]<-Truepositives[z]/(Truepositives[z]+FalsePositives[z])
  F_One_Half[z]<-(1.5*Precision[z]*Sensitivity[z])/(.25*Precision[z]+Sensitivity[z])
  F_One[z]<-(2*Precision[z]*Sensitivity[z])/(Precision[z]+Sensitivity[z])
  F_Two[z]<-(5*Precision[z]*Sensitivity[z])/(4*Precision[z]+Sensitivity[z])
  FalsePositiveRate[z]<-FalsePositives[z]/(TrueNegatives[z]+FalsePositives[z])
}
Characterization<-cbind.data.frame(threshold,Truepositives,Falsenegatives,Sensitivity,TrueNegatives,FalsePositives,One_Minus_Specificity,Accuracy,MCC,EER,Precision,FalsePositiveRate,F_One_Half,F_One,F_Two)

positiveheader<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,"RPMS","PMS")
positivewithscores<-rbind.data.frame(positiveheader,positivewithscores)

negativeheader<-c("Substrate","RPMS","PMS")
colnames(NegativeWithScores)<-negativeheader

write.table(x=c("Characterzation Table"),file = FILENAME2, col.names = FALSE,row.names = FALSE, append = TRUE,sep = ",")
write.table(Characterization,file = FILENAME2, col.names = TRUE,row.names = FALSE, append = TRUE,sep = ",")
#write the output files


options(warn = oldw)
#reset the warning option