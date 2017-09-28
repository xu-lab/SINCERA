# modification of sSeq for Drop-seq and 10x data
# getNormFactor, getT, 

# size factor determined by LibSize/median(LibSize)

getNormFactor.2 <- function(countsTable) {
  
  libsize <- colSums(countsTable)
  sf <- libsize/median(libsize)
  names(sf) <- colnames(countsTable)
  return(sf)
  
}

# replace getNormFactor with getNormFactor.2

getT.2 <- function(
    countsTable, 
    sizeFactors=NULL, 
    q.vec=NULL, 
    plotASD=FALSE, 
    numPart=1, 
    propForSigma=c(0,1), 
    verbose=TRUE, 
    shrinkTarget=NULL, 
    shrinkQuantile=NULL, 
    shrinkVar=FALSE, 
    eSlope=0.05, 
    disp=NULL, 
    dispXX=NULL, 
    normalize=FALSE, 
    lwd1=4.5, 
    cexlab1=1.2
){
    if(!is.null(countsTable)) {counts=as.matrix(countsTable)}
    if(is.null(countsTable) & is.null(disp)){
        stop("at least provide the initial dispersion estimates.")
    }
    if(is.null(sizeFactors) & !is.null(countsTable)) { 
        #mg: sizeFactors=getNormFactor(countsTable)    
		sizeFactors=getNormFactor.2(countsTable) 		
    }
    
    if(is.null(eSlope)){
        eSlope=0.002
    }else{
        if(length(eSlope)>1 & verbose) 
            print("Note: only the first value in eSlope is used for tests.")
    }
		
    allAdjDisp=list()	
    if(is.null(disp) & !is.null(countsTable)){
        normc=as.matrix(t(t(counts)/sizeFactors))
        normc.m=rowMeans(normc)
        normc.v=rowVars(as.matrix(t(t(counts)/sqrt(sizeFactors))))
        s_st=mean(1/sizeFactors)
        disp=(normc.v - normc.m )/(normc.m)^2 
        normc.m[normc.m<=0]=1
        log.normc.m=log(normc.m)
    } else if (!is.null(dispXX)){
        normc.m=dispXX 
        normc.m[normc.m<=0]=1 
        log.normc.m=log(normc.m)
    } else {    
        normc.m=NULL 
        log.normc.m=NULL 
    }
    
    if(shrinkVar&is.null(disp)){
        disp=normc.v 
        if(verbose) 
            print("Shrinkage estimates on variance are used.")
    } else {  
        if(verbose)	
            print("Shrinkage estimates on dispersion are used for the tests.")
    }
    disp[is.na(disp) ]=0
    disp[disp<=0]=0  
    
    if(numPart==1){  
        disp.m= mean(disp)
        asd.mle=round(mean((disp-disp.m)^2, na.rm=T), 4)
        rg.xx=quantile(disp[is.finite(disp)], prob=c(0.05,0.995))
        xx=seq(rg.xx[1], rg.xx[2], length.out=200)
        asd=rep(0, length(xx))    
        for(i in 1:length(xx)){
            allAdjDisp[[i]]=equalSpace(disp, log.normc.m, 1, 
                propForSigma=propForSigma,  shrinkTarget=xx[i], vb=FALSE) 
            allAdjDisp[[i]]=pmax(1e-8, allAdjDisp[[i]])  
            names(allAdjDisp[[i]])=1:length(disp)
            asd[i]=mean((allAdjDisp[[i]]- disp)^2, na.rm=T)
        }
        diff.q=diff.asd=rep(0, length(asd))
        maxASD=max(asd, na.rm=T)
        maxASD.pnt=which( asd == maxASD)
        for ( i in 1:length(asd) ){ 
            diff.asd[i]=maxASD - asd[i]
            diff.q[i]=xx[maxASD.pnt] - xx[i] 
        }
        numAdjPoints=6
        len.asd=length(asd) - numAdjPoints + 1
        slope1=rep(1, len.asd)
        if(normalize){
            xx1=xx/sd(xx)
            yy1=asd/sd(asd)   	
            eSlope=eSlope*5
        } else {
            xx1=xx
            yy1=asd
        }
        for (i in 1:len.asd){ 
            slope1.xx=xx1[ i:(i+numAdjPoints-1) ]
            slope1.yy=yy1[ i:(i+numAdjPoints-1) ]
            slope1[i]=cov(slope1.xx, slope1.yy)/var(slope1.xx)
        }
        maxSlope1=max(abs(slope1))
        maxSlope1.pnt=which(abs(slope1) == maxSlope1)
        sub.slope1=abs(slope1)[maxSlope1.pnt:len.asd]
        sub.diff.asd=diff.asd[maxSlope1.pnt:length(diff.asd)] 
        pred.diff=matrix(NA, nrow=length(sub.diff.asd), ncol=numAdjPoints)
        for(i in 1:length(sub.diff.asd)){
            for(j in 1:numAdjPoints ){ 
                if( i-j >= 0){
                    pred.diff[i,j]=sub.diff.asd[i]/sub.slope1[i-j+1]
                }
            }
        }
        max.pred= max(pred.diff, na.rm=T)
        max.rowInd=which(apply(pred.diff, 1, max, na.rm=T)==max.pred) 
        temp.max.pnt=max.rowInd+maxSlope1.pnt-1-ceiling(numAdjPoints/2)
        max.pnt=rep(0, length(eSlope))
        for(k in 1:length(eSlope)){
            max.pnt[k]=temp.max.pnt[1] 
            while(-slope1[max.pnt[k]-ceiling(numAdjPoints/2)]<eSlope[k]&
                slope1[max.pnt[k]-ceiling(numAdjPoints/2)]<0
            ){
                max.pnt[k]=max.pnt[k] - 1
            }
        }
        if(!is.null(shrinkQuantile)){
            max.pnt=1 
            q.vec=shrinkQuantile
            adjDisp1=equalSpace(disp, log.normc.m, numPart, 
                propForSigma=propForSigma, shrinkQuantile=q.vec[1], vb=FALSE) 
            adjDisp1=pmax(1e-8, adjDisp1)  
            names(adjDisp1)=1:length(disp)
            asd.target=mean((adjDisp1- disp)^2, na.rm=T)   
            target=round(quantile(disp, prob=q.vec[1]),3) 
        }
   
        if(!is.null(shrinkTarget)){
            max.pnt=1 
            disp.tm=c(disp[!is.na(disp)], shrinkTarget[1])
            q.vec=round(rank(disp.tm)[disp.tm==shrinkTarget[1]]/
                length(disp[!is.na(disp)]), 3)
            adjDisp1=equalSpace(disp, log.normc.m, numPart, 
                propForSigma=propForSigma, shrinkQuantile=q.vec[1], vb=FALSE) 
            adjDisp1=pmax(1e-8, adjDisp1)  
            names(adjDisp1)=1:length(disp)
            asd.target=mean((adjDisp1- disp)^2, na.rm=T)
            target=shrinkTarget[1]
        }

        if(is.null(shrinkQuantile) & is.null(shrinkTarget)){
            target=asd.target=q.vec=rep(0,length(eSlope))
            for(k in 1:length(eSlope)){
                target[k]=xx[max.pnt[k]][1]
                asd.target[k]=asd[max.pnt[k]]
                disp.tm=c(disp[!is.na(disp)], target[k])
                q.vec[k]=round(rank(disp.tm)[disp.tm==target[k]]/
                    length(disp[!is.na(disp)]), 3)
            }
        }
   
        if (plotASD) {
            tt1="ASD vs shrinkage target"
            plot(asd~xx, cex=0.5, type="l", lwd=lwd1, main="", ylab='', 
                xlab='', axes=FALSE)
            mtext(text=expression(paste("ASD(",xi, ")", sep='')), side=2, 
                padj=-0.9, cex=cexlab1)
            mtext(text=expression(xi), side=1, padj=1.3, cex=cexlab1)
            axis(1, padj=-1, cex.axis=1.2) 
            axis(2, padj=1, cex.axis=1.2)
            if(!is.null(shrinkTarget) | !is.null(shrinkQuantile)){
                abline(h=asd.target, lty=2, col="blue")
                abline(v=target, lty=2, col="blue")     	
                mtext(text=round(asd.target, 3), side=4, at=asd.target, las=2)
                mtext(text=round(target,3), side=3, at=target, las=0)
            } else {
                prev.max.pnt=0
                for(k in 1:length(eSlope)){
                    if(prev.max.pnt==max.pnt[k]) 
                        next       	
                    abline(h=asd.target[k], lty=2, col="gray20")
                    abline(v=target[k], lty=2, col="gray20")
                    mtext(text=round(asd.target[k], 3), side=4, 
                        at=asd.target[k], las=2)
                    if(k==1){
                        leg.k=as.expression( 
                            bquote(hat(italic(xi))~.(paste("=", 
                            round(target[k],3),sep='')) ))
                    } else {
                        leg.k=target[k]
                    }
                    mtext(text=leg.k, side=3, at=round(target[k],3), las=0)
                    prev.max.pnt=max.pnt[k]
                }
            }
            abline(h=asd.mle, lty=1, col="gray40")
            mtext(text=asd.mle, side=4, at=asd.mle, las=2)
        }
   
        if(verbose){
            if(!is.null(shrinkQuantile) | !is.null(shrinkTarget) ){
                print(paste("The selected shrink target is", target[1]))
                print(paste("The selected shrink quantile is", q.vec[1]))
            } else { 
                print(paste("The shrink target is", target[1])) 
                print(paste("The shrink quantile is", q.vec[1])) 
            }
        }
        return(list(q=q.vec[1], target=target[1]))
   
    } else if (numPart>1) {
        if(is.null(log.normc.m)){
            stop("Error in getT: log.normc.m can not be NULL.")
        }
        out=getTgroup(y=disp, x=log.normc.m, numPart=numPart, plotASD=plotASD,
            verbose=verbose, eSlope=eSlope, lwd1=lwd1, cexlab1=cexlab1)
        return(out)	
    } else {
        stop("Error: numPart must be a non-negative integer.")
    }
}

# replace getT and getNormFactor with getT.2 and getNormFactor.2, respectively

nbinomTestForMatricesSH.2 <- function(
    countsA, 
    countsB, 
    sizeFactorsA, 
    sizeFactorsB, 
    numPart=1, 
    SHonly=FALSE, 
    propForSigma=c(0,1), 
    shrinkTarget=NULL, 
    shrinkQuantile=NULL, 
    cLA, 
    cLB, 
    contrast=NULL, 
    keepLevelsConsistant=TRUE,
    useMMdisp=FALSE,
    shrinkVariance=FALSE,
    pairedDesign=FALSE, 
    pairedDesign.dispMethod="per-pair", 
    useFisher=FALSE,
    Dispersions=NULL,
    eSlope=NULL,
    plotASD=FALSE,
    lwd_ASD=4.5, 
    cex_ASD=1.2
){
    cl.nm=sort(unique(c(cLA, cLB))) 
    cntA=as.matrix(countsA) 
    cntB=as.matrix(countsB) 
    sfA=sizeFactorsA 
    sfB=sizeFactorsB 
    s.m=1/mean(1/c(sfA, sfB))  
    s_st=mean(1/c(sfA, sfB))    
    s.tilde1=sum(c(sfA, sfB)) 
    normA=t(t(cntA)/sfA) 
    normB=t(t(cntB)/sfB)
    norm=cbind(normA, normB)
    norm.m=rowMeans(norm)
    norm.v=rowVars(norm)
    norm.mA=rowMeans(normA)
    norm.vA=rowVars(normA)
    norm.mB=rowMeans(normB) 
    norm.vB=rowVars(normB)    
    orig.shrinkQuantile=shrinkQuantile
    orig.shrinkTarget=shrinkTarget
    adj.getT=FALSE
    ll=list()    
    for(L in 1:length(cl.nm)) {
        countsA=as.matrix(cntA[,cLA==cl.nm[L]]) 
        nA=ncol(countsA)
        countsB=as.matrix(cntB[,cLB==cl.nm[L]]) 
        nB=ncol(countsB)
        if( is.null(orig.shrinkQuantile) & is.null(shrinkTarget) & 
            is.null(Dispersions)
        ){
            if(length(cl.nm)>1) 
                print(paste("Get shrinkage target at level", cl.nm[L]))
            countsTable1=cbind(countsA, countsB)
            #mg: sf1=getNormFactor(countsTable1)
			sf1=getNormFactor.2(countsTable1)
            #mg: out.getT=getT(countsTable1, sf1, numPart=numPart, plotASD=FALSE, propForSigma=propForSigma, verbose=FALSE, shrinkVar=shrinkVariance, eSlope=eSlope, disp=NULL, wd1=lwd_ASD, cexlab1=cex_ASD)
			out.getT=getT.2(countsTable1, sf1, numPart=numPart, plotASD=FALSE, 
                propForSigma=propForSigma, verbose=FALSE, 
                shrinkVar=shrinkVariance, eSlope=eSlope, disp=NULL, 
                lwd1=lwd_ASD, cexlab1=cex_ASD) 
            shrinkTarget=out.getT$target
            if(min(out.getT$target)>=0.8 & min(out.getT$q)>=0.965){
                adj.getT=TRUE
            }
        }
        N=nA + nB 
        rA=nA/N  
        rB=nB/N
        cnt=cbind(countsA, countsB) 
        ng=nrow(cnt)
        kAs=rowSums(cbind(countsA)) 
        kBs=rowSums(cbind(countsB))
        conds0=c(rep('A',nA),rep('B',nB))
        countsTable0=cbind(countsA, countsB)
        #mg: sizeFactors=getNormFactor(countsTable0)  
		sizeFactors=getNormFactor.2(countsTable0)    		
        sizeFactorsA=sizeFactors[conds0=="A"]
        sizeFactorsB=sizeFactors[conds0=="B"]
        normcA=t(t(countsA)/sizeFactorsA)
        normcB=t(t(countsB)/sizeFactorsB)   
        normc=as.matrix(cbind(normcA, normcB), ncol=nA+nB)
        normc.m=rowMeans(normc)
        normc.v=rowVars(normc)
        s=sum(sizeFactors)
        s.mL=1/mean(1/sizeFactors)
        s.tilde=sum(sizeFactors)
        dispc=s.mL*(normc.v - normc.m*s_st )/(normc.m)^2
        dispc[is.na(dispc) ]=0
        dispc[dispc<=0]=0 
        xx=normc.m
        xx[xx<=0]=1
        if(shrinkVariance){
            normc.v[is.na(normc.v)]=0
            normc.v=equalSpace(normc.v, log(xx), numPart, 
                propForSigma=propForSigma, shrinkTarget=shrinkTarget, 
                shrinkQuantile=shrinkQuantile, vb=FALSE)
            adjdispc=s.m*(normc.v-normc.m*s_st )/(normc.m)^2
            adjdispc[is.na(adjdispc) ]=0 
            adjdispc[adjdispc<=0]=0
        } else {
            if(is.null(Dispersions)){
                adjdispc=equalSpace(dispc, log(xx), numPart, 
                    propForSigma=propForSigma, shrinkTarget=shrinkTarget, 
                    shrinkQuantile=shrinkQuantile, vb=FALSE) 
            } else {
                adjdispc=Dispersions
            }
        }
        sumDispsc=adjdispc 
        names(sumDispsc)=rownames(countsA)
        ll[[L]]=list(countsA=countsA, countsB=countsB, 
                    normcA=normcA, normcB=normcB,
                    sumDisps=adjdispc, disp=dispc, sfA=sizeFactorsA, 
                    sfB=sizeFactorsB, mus=normc.m, s=s, 
                    sA=sum(sizeFactorsA), sB=sum(sizeFactorsB),
                    s.tilde=s.tilde, nA=nA, nB=nB, s.m=s.mL, 
                    counts=cnt, kAs=kAs, kBs=kBs,  rA=rA, rB=rB)
    }
    s.m1=1
    for(L in 1:length(cl.nm)){ 
        s.m1=max(ll[[L]]$s.m, s.m1) 
    }
    shrinkQuantile=orig.shrinkQuantile
    shrinkTarget=orig.shrinkTarget    
    verbose=TRUE
    if(pairedDesign){
        if(pairedDesign.dispMethod=="per-pair"){
            print(paste("For paired design, the per-pair dispersion", 
                "estimates are used and shrunk separately."))
            verbose=FALSE
            disp=sumDisps=NULL
            for(L in 1:length(cl.nm)){
                disp=cbind(disp, ll[[L]]$disp)
                sumDisps=cbind(sumDisps, ll[[L]]$sumDisps)
            }
        } else if (pairedDesign.dispMethod=="pooled"){
            print(paste("For paired design, the aveaged dispersion", 
                "estimates across paires are used."))
            disp.pair=NULL
            for(L in 1:length(cl.nm)){
                disp.pair=cbind(disp.pair, ll[[L]]$disp)
            }
            disp=rowMeans(disp.pair, na.rm=T)
        } else {
            stop(paste("Error in defining the method of dispersion", 
                "estimates for paired design. \nPlease input 'per-pair'", 
                "or 'pooled'."))
        }
    } else {
        disp=s.m*(norm.v-norm.m*s_st)/(norm.m)^2
        disp[ is.na(disp) ]=0  
        disp[disp<=0]=0 
    }
    if(shrinkVariance){
        xx=norm.m 
        xx[xx<=0]=1    
        norm.v[is.na(norm.v)]=0  
        if(is.null(shrinkQuantile) & is.null(shrinkTarget)){
            #mg: shrinkTarget=getT(countsTable=NULL, sizeFactors=NULL, numPart=numPart, propForSigma=propForSigma, verbose=FALSE,plotASD=plotASD, shrinkVar=shrinkVariance, eSlope=eSlope, disp=norm.v, dispXX=norm.m, lwd1=lwd_ASD, cexlab1=cex_ASD)$target
			shrinkTarget=getT.2(countsTable=NULL, sizeFactors=NULL, numPart=numPart, propForSigma=propForSigma, verbose=FALSE,plotASD=plotASD, shrinkVar=shrinkVariance, eSlope=eSlope, disp=norm.v, dispXX=norm.m, lwd1=lwd_ASD, cexlab1=cex_ASD)$target
        }
        norm.v=equalSpace(norm.v, log(xx), numPart, 
            propForSigma=propForSigma, shrinkTarget=shrinkTarget, 
            shrinkQuantile=shrinkQuantile, vb=FALSE)    
        adjdisp=s.m*( norm.v - norm.m*s_st )/(norm.m)^2
        adjdisp[ is.na(adjdisp) ]=0  
        adjdisp[adjdisp<=0]=0
        sumDisps=adjdisp
    } else {
        if(useMMdisp){
            sumDisps=disp
            adjdisp=rep(NA, length(sumDisps))
        } else {  
            if(is.null(shrinkQuantile) & is.null(shrinkTarget) & 
                is.null(Dispersions)
            ){
                #mg: out.getT=getT(countsTable=NULL, sizeFactors=NULL, numPart=numPart, plotASD=plotASD, propForSigma=propForSigma, verbose=verbose, shrinkVar=shrinkVariance, eSlope=eSlope, disp=disp, dispXX=norm.m, lwd1=lwd_ASD, cexlab1=cex_ASD)
				out.getT=getT.2(countsTable=NULL, sizeFactors=NULL, numPart=numPart, plotASD=plotASD, propForSigma=propForSigma, verbose=verbose, shrinkVar=shrinkVariance, eSlope=eSlope, disp=disp, dispXX=norm.m, lwd1=lwd_ASD, cexlab1=cex_ASD)
                shrinkTarget=out.getT$target
                if(min(out.getT$target)>=0.8 & min(out.getT$q)>=0.965){
                    adj.getT=TRUE
                }
            }
            if(!pairedDesign | 
                (pairedDesign & pairedDesign.dispMethod=="pooled")
            ){ 
                xx=norm.m
                xx[xx<=0]=1      
                if(is.null(Dispersions)){
                    adjdisp=equalSpace(disp, log(xx), numPart, 
                        propForSigma=propForSigma, vb=FALSE, 
                        shrinkTarget=shrinkTarget, 
                        shrinkQuantile=shrinkQuantile)     
                } else {
                    adjdisp=Dispersions
                }
                sumDisps=adjdisp
            }
        }
    }
    if(length(ll)>1){ 
        temCount=ll[[1]]$mus
        for(L in 2:length(ll)){
            temCount=cbind(temCount, ll[[L]]$mus)
        }
        lm=rowMeans(log(temCount))
        int.func=function(x1,lm){ 
            exp(median((log(x1)-lm)[is.finite(lm)]))
        }
        s_L=apply(temCount, 2, int.func, lm)
        s_L1 =sum(s_L)
        s_L.m=s_L1/length(s_L)
    } else {
        s_L=1 
        s_L.m=1
    }
    if(SHonly) {
        return(data.frame(SH=sumDisps, raw=disp, mus=norm.m))
    }
    if(!is.null(Dispersions)){
        if(length(Dispersions)!=length(norm.m)){
            stop(paste("Please let the length of input Dispersion", 
                "equal to the number of rows in countsTable, or", 
                "set it as NULL."))
        }
        print("The known dispersion values are used.")
        sumDisps=Dispersions
    }    
    perc=c(1,3,5,7,9, 10)
    progress=round(perc/10 * ng, 0) 
    progress.id=1
    pval=rep(1, ng)
    if(pairedDesign){
        for(i in 1:ng){
            if(progress.id<length(progress)&i==progress[progress.id]){
                progress.id=progress.id + 1
                print(paste(perc[progress.id]*10, "% processed.", sep=''))
            }
            log.p=0
            for(L in 1:length(cl.nm)){
                kA1=ll[[L]]$kAs[i] 
                kB1=ll[[L]]$kBs[i] 
                norm.mA1=s_L.m*norm.m[i]*ll[[L]]$sA 
                norm.mB1=s_L.m*norm.m[i]*ll[[L]]$sB 
                sA1=s_L.m*ll[[L]]$s.tilde
                sB1=s_L.m*ll[[L]]$s.tilde
                ks=kA1 + kB1       
                if(pairedDesign.dispMethod=="per-pair"){
                    allps1=dnbinom( 0:ks, mu= norm.mA1,  
                        size=sA1/ll[[L]]$sumDisps[i])*dnbinom(ks:0, 
                        mu=norm.mB1, size=sB1/ll[[L]]$sumDisps[i] ) 
                    pobs1=dnbinom(kA1, mu=norm.mA1, 
                        size=sA1/ll[[L]]$sumDisps[i])*dnbinom(kB1,
                        mu=norm.mB1, size=sB1/ll[[L]]$sumDisps[i] )      
                } else {
                    allps1=dnbinom( 0:ks, mu= norm.mA1,  
                        size=sA1/sumDisps[i])*dnbinom(ks:0, 
                        mu=norm.mB1, size=sB1/sumDisps[i] ) 
                    pobs1=dnbinom(kA1, mu=norm.mA1, 
                        size=sA1/sumDisps[i])*dnbinom(kB1, 
                        mu=norm.mB1, size=sB1/sumDisps[i] )      
                }
                sumsel=sum(allps1[allps1<=pobs1], na.rm=T)
                sumall=sum(allps1, na.rm=T) 
                log.p1=log(min(1, sumsel/sumall))
                if(is.finite(log.p1)) {
                    log.p=log.p + log.p1
                } 
            }
            if(useFisher){
                chi1=-2 * log.p
                pval[i]=1-pchisq(chi1, length(cl.nm)*2)        
            } else {
                pval[i]=exp(log.p)
            }
        } 
        pval[is.na(pval)]=1 
        return(list(pval=pval, dispMM=disp, dispSH=sumDisps, mu=norm.m))        
    }  
    
    if(is.null(contrast)){    
        for(i in 1:ng){
            if(progress.id<length(progress) & 
                i==progress[progress.id] 
            ){
                progress.id=progress.id + 1
                print(paste(perc[progress.id]*10,"% processed.", sep=''))
            }        
            sumsel=sumall=0  
            sign.diff=0 
            for(L in 1:length(cl.nm)){
                ks=(ll[[L]]$kAs[i] + ll[[L]]$kBs[i])
                nks=length(0:ks) 
                s=ll[[L]]$s 
                rA=ll[[L]]$rA 
                rB=ll[[L]]$rB     
                if(ll[[L]]$nA>1 & ll[[L]]$nB>1){    
                    norm.mA1=s_L.m*norm.m[i]*ll[[L]]$sA
                    norm.mB1=s_L.m*norm.m[i]*ll[[L]]$sB
                    sA1=s_L.m*ll[[L]]$sA*2
                    sB1=s_L.m*ll[[L]]$sB*2
                } else {
                    norm.mA1=s_L.m*norm.m[i]*ll[[L]]$s*ll[[L]]$rA
                    norm.mB1=s_L.m*norm.m[i]*ll[[L]]$s*ll[[L]]$rB              
                    sA1=s_L.m*ll[[L]]$s.tilde
                    sB1=s_L.m*ll[[L]]$s.tilde
                }
                allps1=dnbinom(0:ks, mu=norm.mA1, size=sA1/sumDisps[i])*
                    dnbinom(ks:0, mu=norm.mB1, size=sB1/sumDisps[i] )  
                pobs1=dnbinom(ll[[L]]$kAs[i], mu=norm.mA1, 
                    size=sA1/sumDisps[i])*dnbinom(ll[[L]]$kBs[i], 
                    mu=norm.mB1, size=sB1/sumDisps[i]) 
                sumsel=sumsel+ sum(allps1[allps1<=pobs1], na.rm=T)
                sumall=sumall+ sum(allps1, na.rm=T)
                sign.diff=sign.diff + 
                    (mean(ll[[L]]$normcA[i])<mean(ll[[L]]$normcB[i]))
            }
            if(keepLevelsConsistant & sign.diff>0&
                sign.diff<length(cl.nm)
            ){
                pval[i]=1
            } else {
                pval[i]=min(1, sumsel/sumall)        
            }
        }
        pval[is.na(pval)]=1
        if(adj.getT){
            rk=rank(pval, ties.method="min")
            wt=rk/max(rk)
            pval=pval*wt
        }       
    } else { 
        for(i in 1:ng){
            if (progress.id<length(progress)&i==progress[progress.id]){
                progress.id=progress.id + 1
                print(paste(perc[progress.id]*10,"% processed.", sep=''))
            }
            norm.mA1 =norm.mB1=kA1=kB1=sA1.sz=sB1.sz=
                sA1=sB1=sumDisps1= 0
            for(L in 1:length(cl.nm)){
                kA1=ll[[L]]$kAs[i] + kA1
                kB1=ll[[L]]$kBs[i] + kB1
                norm.mA1=ll[[L]]$mus[i]*ll[[L]]$sA+norm.mA1 
                norm.mB1=ll[[L]]$mus[i]*ll[[L]]$sB+norm.mB1 
                sA1=s_L.m*ll[[L]]$s.tilde + sA1
                sB1=s_L.m*ll[[L]]$s.tilde + sB1
            }
            sumDisps1=sumDisps[i] 
            ks= kA1 + kB1
            nks=length(0:ks)         
            allps1=dnbinom(0:ks, mu=norm.mA1, size=sA1/sumDisps1)*
                dnbinom(ks:0, mu=norm.mB1, size=sB1/sumDisps1)   
            pobs1=dnbinom(kA1, mu=norm.mA1, size=sA1/sumDisps1)*
                dnbinom(kB1, mu=norm.mB1, size=sB1/sumDisps1)  
            sumsel=sum(allps1[allps1<=pobs1], na.rm=T)
            sumall=sum(allps1, na.rm=T)
            pval[i]=min(1, sumsel/sumall)
        }
    }
    pval[is.na(pval)]=1    
    return(list(pval=pval, dispMM=disp, dispSH=sumDisps, mu=norm.m))
}


# replace getNormFactor and nbinomTestForMatricesSH with getNormFactor.2 and nbinomTestForMatricesSH.2, respectively

nbTestSH.2 <- function( 
    countsTable, 
    conds,
    condA="A", 
    condB="B",
    numPart=1, 
    SHonly=FALSE, 
    propForSigma=c(0, 1), 
    shrinkTarget=NULL, 
    shrinkQuantile=NULL, 
    plotASD=FALSE, 
    coLevels=NULL, 
    contrast=NULL, 
    keepLevelsConsistant=FALSE, 
    useMMdisp=FALSE, 
    addRawData=FALSE, 
    shrinkVariance=FALSE,
    pairedDesign=FALSE,
    pairedDesign.dispMethod="per-pair", 
    useFisher=FALSE,
    Dispersions=NULL,
    eSlope=0.05,
    lwd_ASD=4.5, 
    cex_ASD=1.2
){
    if(is.null(coLevels)){
        coLevels=data.frame(exp=rep(1,ncol(countsTable)))
    }
    for(j in 1:ncol(coLevels)){
        cL=apply(coLevels, 1, paste, collapse="_")
    }
    cL.tb=table(cL)       
    if( sum(cL.tb<2)>0 ){
       stop("Errors in 'coLevels'. All the levels must have paired comparisons.")
    }
    if(!is.null(contrast)){
        if(length(contrast)!=ncol(countsTable)){
            stop(paste("Error: the length of contrast vector must equal", 
                "to the number of columns in countsTable."))
        }
        if(length(unique(conds[contrast>0]))!=1){
            stop(paste("Error: this package is currently only available", 
                "for the contrast between conditions, not the contrast", 
                "within conditions. Please revise the contrast vector,", 
                "such as c(1,1,-1,-1) for cond=c('A','A', 'B, 'B'),", 
                "instead of c(1,-1,1,-1)."))
        }
        countsTable=countsTable[, contrast!=0]
        conds=conds[contrast!=0]
        cL=cL[contrast!=0]  
        contrast= contrast[contrast!=0]
        contrast.mat= t(matrix(contrast, nrow=length(contrast), 
            ncol=nrow(countsTable)))
        if(sum(abs(contrast)!=1)>0){
            countsTable=round(countsTable * abs(contrast.mat))
        }
    }
    #mg: sf=getNormFactor(countsTable)
	sf=getNormFactor.2(countsTable)

    colA=conds == condA
    colB=conds == condB

    if(sum(colB[1:sum(colA)])>0){
        stop(paste("re-order the columns of the countsTable so that all", 
            "samples in the same condition are in the adjacent columns.", 
            "For example, 'A A B B' is good, but 'A B A B' is bad."))
    }
    
    counts=as.matrix(countsTable)
    ng=nrow(counts)
    cntA=matrix(counts[, colA], ncol=sum(colA)) 
    cntB=matrix(counts[, colB], ncol=sum(colB))
        
    if(SHonly){
        #mg: disp=nbinomTestForMatricesSH(countsA=cntA, countsB=cntB, sizeFactorsA=sf[colA], sizeFactorsB=sf[colB],numPart=numPart, SHonly=TRUE, propForSigma= propForSigma, shrinkTarget=shrinkTarget, shrinkQuantile=shrinkQuantile,cLA=cL[colA], cLB=cL[colB], keepLevelsConsistant, useMMdisp,shrinkVariance=shrinkVariance, eSlope=eSlope, plotASD=plotASD, lwd_ASD=lwd_ASD, cex_ASD=cex_ASD)
			
		disp=nbinomTestForMatricesSH.2(countsA=cntA, countsB=cntB, 
            sizeFactorsA=sf[colA], sizeFactorsB=sf[colB],
            numPart=numPart, SHonly=TRUE, propForSigma= propForSigma, 
            shrinkTarget=shrinkTarget, shrinkQuantile=shrinkQuantile,
            cLA=cL[colA], cLB=cL[colB], keepLevelsConsistant, useMMdisp,
            shrinkVariance=shrinkVariance, eSlope=eSlope, plotASD=plotASD, 
            lwd_ASD=lwd_ASD, cex_ASD=cex_ASD)
        return(disp)
    } else {
        t1=Sys.time()
        #mg: pval0=nbinomTestForMatricesSH(countsA=cntA, countsB=cntB, sizeFactorsA=sf[colA], sizeFactorsB=sf[colB],numPart=numPart, SHonly=FALSE, propForSigma= propForSigma, shrinkTarget=shrinkTarget, shrinkQuantile=shrinkQuantile, cLA=cL[colA], cLB=cL[colB], contrast=contrast,  keepLevelsConsistant=keepLevelsConsistant, useMMdisp,shrinkVariance=shrinkVariance, pairedDesign=pairedDesign, pairedDesign.dispMethod=pairedDesign.dispMethod, useFisher=useFisher, Dispersions=Dispersions, eSlope=eSlope, plotASD=plotASD,  lwd_ASD=lwd_ASD, cex_ASD=cex_ASD)
		
		pval0=nbinomTestForMatricesSH.2(countsA=cntA, countsB=cntB, 
            sizeFactorsA=sf[colA], sizeFactorsB=sf[colB],
            numPart=numPart, SHonly=FALSE, propForSigma= propForSigma, 
            shrinkTarget=shrinkTarget, shrinkQuantile=shrinkQuantile, 
            cLA=cL[colA], cLB=cL[colB], contrast=contrast,  
            keepLevelsConsistant=keepLevelsConsistant, useMMdisp,
            shrinkVariance=shrinkVariance, pairedDesign=pairedDesign, 
            pairedDesign.dispMethod=pairedDesign.dispMethod, 
            useFisher=useFisher, Dispersions=Dispersions, 
            eSlope=eSlope, plotASD=plotASD,  lwd_ASD=lwd_ASD, 
            cex_ASD=cex_ASD)
        print(Sys.time()-t1)
        dispMM=pval0$dispMM
        dispSH=pval0$dispSH
        Mean=pval0$mu
        pval=pval0$pval
    }
    cl.nm=sort(unique(c(cL[colA], cL[colB])))
    rM.A=rowMeans(as.matrix(counts[,colA])) 
    rM.B=rowMeans(as.matrix(counts[,colB]))
    l2f=log2(rM.A/rM.B)
    
    rs=data.frame(  
        Mean=Mean,
        rawMeanA=rM.A,
        rawMeanB=rM.B, 
        rawLog2FoldChange=l2f, 
        dispMM=dispMM,
        dispSH=dispSH,
        pval=pval, 
        stringsAsFactors=FALSE)
    rownames(rs)=rownames(countsTable)
    if(addRawData){
        rs=data.frame(rs, countsTable)
    }
    return(rs)
}

