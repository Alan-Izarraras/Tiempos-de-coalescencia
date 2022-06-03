UpperThreshold = 60
Limit = UpperThreshold + 2.5
Runs = 60/5
RunsPlusTwo = Runs + 3

#### These two lists define the grid of the alpha and gamma parameters
AlphaGrid <- 0.02*1:10
GammaGrid <- 2*1:10

TwoNsValues <- 0:Runs * 5 + 2.5

Table <- matrix(ncol=RunsPlusTwo,nrow=0)

for (j in AlphaGrid){
    for (k in GammaGrid){
		Probability <- 0
		Row <- c(j,k)
		for (i in TwoNsValues){
			
			if ( i == 2.5){
				Probability <- pgamma(2.5,j,1/k)
				Row <- c(Row,Probability)
			}else if (i==Limit){
				Probability <- ( 1 - pgamma(i-5,j,1/k) )
				Row <- c(Row,Probability)
			}else{
				Probability <-(pgamma(i,j,1/k) - pgamma(i-5,j,1/k))
				Row <- c(Row,Probability)
			}
			
		}
		Table <- rbind(Table,Row)
#print (Probability)
    }
}

colnames(Table) <- c("alpha", "gamma", "P(2Ns < 2.5)", "P(2.5 < 2Ns < 7.5)", "P(7.5 < 2Ns < 12.5)", "P(12.5 < 2Ns < 17.5)", "P(17.5 < 2Ns < 22.5)", "P(22.5 < 2Ns < 27.5)", "P(27.5 < 2Ns < 32.5)", "P(32.5 < 2Ns < 37.5)", "P(37.5 < 2Ns < 42.5)", "P(42.5 < 2Ns < 47.5)", "P(47.5 < 2Ns < 52.5)", "P(52.5 < 2Ns < 57.5)", "P( 2Ns > 57.5)")

write.table(Table,file="TableOfProbabilities.txt",row.names=FALSE,col.names=FALSE,sep="\t")

### Show that the probabilities are equal to 1
for (k in 1:100){
    print(sum(Table[k,3:15]))
}

