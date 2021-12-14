# code from Harmon & Glor 2010

# Calculates probabilities for phylogenetic permutations; from Lapointe and Garland 2001
phyloProb<-function(phy, k=2) {
	pd<-cophenetic.phylo(phy)
	pdr<-pd/max(pd)
	s<-k-pdr
	p<-s/rowSums(s)
	p
	}

# Permutes species according to phylogentic tree; returns tip names in permuted order
phyloPermute<-function(phy, k=1) {
	p<-phyloProb(phy, k)
	tt<-rownames(p)
	nsp<-length(tt)
		order<-sample(1:nsp, replace=F)
		ttnew<-character(nsp)
		cpm<-p
		for(j in order[-nsp]) {
			cpm<-cpm/rowSums(cpm)
			rr<-which(rownames(cpm)==tt[j])
			pp<-cpm[rr,]
			s2<-sample(names(pp), size=1, replace=T, prob=pp)
			slot<-which(tt==s2)
			rc<-which(colnames(cpm)==s2)
			
			ttnew[slot]<-tt[j]
			cpm<-cpm[-rr,-rc]
			}
		ttnew[which(ttnew=="")]<-tt[order[nsp]]

	ttnew
}	

# Mantel test with phylogenetic permutations
phyloMantel <- function(phy,m1,m2,m3=NULL,nperm=999,graph=F,rpt=F) 
{ 
	if(is.null(m3))
	{
		type<-"Two-way Mantel test"
		rm1<-m1
		rm2<-m2
	}
	else
	{
		type<-"Three-way Mantel test"
		rm1 <- mant.resid(m1, m3)
		rm2 <- mant.resid(m2, m3)
		rownames(rm1)<-rownames(m1)
		rownames(rm2)<-rownames(m2)
		colnames(rm1)<-colnames(m1)
		colnames(rm2)<-colnames(m2)
	}
		
	n <- dim(rm1)[1] 
	realz <- mant.zstat(rm1,rm2)
	corr <- mant.r(rm1, rm2) 
	enum<-F 

	nullstats <- numeric(nperm) 
	for (i in (1:nperm)) { 

			nullstats[i] <- mant.zstat(rm1,permPhylo(phy, rm2,n)) 

	} 
	pval <- sum(nullstats>=realz)/(nperm+1);
	if (graph) { 
			plot(density(nullstats),type="l", 
			main="Distribution of Mantel z-statistic", 
			xlab="Z statistic",ylab="# of permutations", 
			sub=paste("Actual z-stat =",round(realz,3), 
			": p<",round(pval,4),",",nperm," permutations")) 
			abline(v=realz) 
	} 
	list(test=type, z.stat=realz*2,p=pval,r = corr) 
	
}


# permute rows and columns of a matrix 
# m1 is a (square) nxn matrix 
# if p is numeric then the pth permutation will be taken, otherwise 
# a random permutation will be taken 
permPhylo <- function(phy, m1,n) 
{ 
s <- phyloPermute(phy) 
trans<-match(s, rownames(m1))
m1[s,s] 
} 