## ----echo=F------------------------------------------------------------------------
options("width"=85)
library(gRim)
ps.options(family="serif")

## ----------------------------------------------------------------------------------
args(dmod)
args(cmod)
args(mmod)

## ----------------------------------------------------------------------------------
data(reinis)
str(reinis)

## ----print=F-----------------------------------------------------------------------
data(reinis)
dm1 <- dmod(list(c("smoke", "systol"), c("smoke", "mental", "phys")), data=reinis)
dm1 <- dmod(~smoke:systol + smoke:mental:phys, data=reinis)
dm1 

## ----------------------------------------------------------------------------------
formula(dm1)
terms(dm1)

## ----print=F-----------------------------------------------------------------------
dm2 <- dmod(~.^2, margin=c("smo","men","phy","sys"),
            data=reinis)
formula(dm2)

## ----print=F-----------------------------------------------------------------------
dm3 <- dmod(list(c("smoke", "systol"), c("smoke", "mental", "phys")),
            data=reinis, interactions=2)
formula(dm3)

## ----fig=T-------------------------------------------------------------------------
plot(dm1)

## ----------------------------------------------------------------------------------
data(carcass)
cm1 <- cmod(~Fat11:Fat12:Fat13, data=carcass)
cm1 <- cmod(~Fat11:Fat12 + Fat12:Fat13 + Fat11:Fat13, data=carcass)
cm1

## ----fig=T-------------------------------------------------------------------------
plot(cm1)

## ----------------------------------------------------------------------------------
data(milkcomp1)
mm1 <- mmod(~.^., data=milkcomp1)
mm1

## ----fig=T-------------------------------------------------------------------------
plot(mm1) ## FIXME: should use different colours for disc and cont variables.

## ----------------------------------------------------------------------------------
###  Set a marginal saturated model:
ms <- dmod(~.^., marginal=c("phys","mental","systol","family"), data=reinis)
formula(ms)
###   Delete one edge:
ms1 <- update(ms, list(dedge=~phys:mental))
formula(ms1)
###   Delete two edges:
ms2<- update(ms, list(dedge=~phys:mental+systol:family))
formula(ms2)
###   Delete all edges in a set:
ms3 <- update(ms, list(dedge=~phys:mental:systol))
formula(ms3)
### Delete an interaction term
ms4 <- update(ms, list(dterm=~phys:mental:systol) )
formula(ms4)

## ----------------------------------------------------------------------------------
###  Set a marginal independence model:
m0 <- dmod(~.^1, marginal=c("phys","mental","systol","family"), data=reinis)
formula(m0)

###  Add three interaction terms:
ms5 <- update(m0, list(aterm=~phys:mental+phys:systol+mental:systol) )
formula(ms5)
###  Add two edges:
ms6 <- update(m0, list(aedge=~phys:mental+systol:family))
formula(ms6)

## ----print=T-----------------------------------------------------------------------
cit <- ciTest(reinis, set=c("systol", "smoke", "family", "phys"))
cit 

## ----------------------------------------------------------------------------------
cit$slice

## ----------------------------------------------------------------------------------
ciTest(reinis, set=c("systol","smoke","family","phys"), method='MC')

## ----print=T-----------------------------------------------------------------------
dm5 <- dmod(~ment:phys:systol + ment:systol:family + phys:systol:smoke,
            data=reinis)

## ----fundamentalfig1,fig.cap="Model for reinis data.", echo=F----------------------
plot(dm5)

## ----------------------------------------------------------------------------------
testdelete(dm5, ~smoke:systol)
testdelete(dm5, ~family:systol)

## ----------------------------------------------------------------------------------
testadd(dm5, ~smoke:mental)

## ----print=T-----------------------------------------------------------------------
ed.in <- getInEdges(ugList(terms(dm5)), type="decomposable")

## ----print=T-----------------------------------------------------------------------
ed.out <- getOutEdges(ugList(terms(dm5)), type="decomposable")

## ----------------------------------------------------------------------------------
args(testInEdges)
args(testOutEdges)

## ----------------------------------------------------------------------------------
testInEdges(dm5, getInEdges(ugList(terms(dm5)), type="decomposable"),
             k=log(sum(reinis)))

## ----fig=T-------------------------------------------------------------------------
dm.sat <- dmod(~.^., data=reinis)
dm.back <- backward(dm.sat)
plot(dm.back)

## ----------------------------------------------------------------------------------
cm.sat <- cmod(~.^., data=carcassall[,1:15])
cm.back <- backward(cm.sat, k=log(nrow(carcass)), type="unrestricted")
plot(cm.back)

## ----fig=T-------------------------------------------------------------------------
dm.i   <- dmod(~.^1, data=reinis)
dm.forw <- forward(dm.i)
plot(dm.forw)

## ----------------------------------------------------------------------------------
dm.s2<-stepwise(dm.sat, details=1)

## ----------------------------------------------------------------------------------
dm.i2<-stepwise(dm.i, direction="forward", details=1)

## ----stepwise01, fig=T, include=F--------------------------------------------------
par(mfrow=c(1,2))
dm.s2
dm.i2
plot(dm.s2)
plot(dm.i2)

## ----------------------------------------------------------------------------------
fix <- list(c("smoke","phys","systol"), c("systol","protein"))
fix <- do.call(rbind, unlist(lapply(fix, names2pairs),recursive=FALSE))
fix
dm.s3 <- backward(dm.sat, fixin=fix, details=1)

## ----------------------------------------------------------------------------------
dm.i3 <- forward(dm.i, fixout=fix, details=1)

## ----stepwise02, fig=T, include=F--------------------------------------------------
par(mfrow=c(1,2))
dm.s3
dm.i3
plot(dm.s3)
plot(dm.i3)

## ----fig=T-------------------------------------------------------------------------
data(mildew)
dm1 <- dmod(~.^., data=mildew)
dm1
dm2 <- stepwise(dm1)
dm2
plot(dm2)

## ----------------------------------------------------------------------------------
ff <- ~la10:locc:mp58:c365+mp58:c365:p53a:a367
mm <- dmod(ff, data=mildew)
plot(mm)

## ----------------------------------------------------------------------------------
dim_loglin(terms(mm), mildew)
dim_loglin_decomp(terms(mm), mildew)

## ----------------------------------------------------------------------------------
data(reinis)
ff    <- ~smoke:mental+mental:phys+phys:systol+systol:smoke
dmod(ff, data=reinis)

## ----------------------------------------------------------------------------------
glist <- rhsFormula2list(ff)
glist
fv1 <- loglin(reinis, glist, print=FALSE)
fv1[1:3]

## ----------------------------------------------------------------------------------
fv2 <- effloglin(reinis, glist, print=FALSE)
fv2[c('logL','nparm','df')]

## ----------------------------------------------------------------------------------
stab <- lapply(glist, function(gg) tableMargin(reinis, gg))
fv3 <- effloglin(stab, glist, print=FALSE)

## ----------------------------------------------------------------------------------
m1 <- loglin(reinis, glist, print=F, fit=T)
f1 <- m1$fit
m3 <- effloglin(stab, glist, print=F, fit=T)
f3 <- m3$fit
max(abs(f1 %a-% f3))

## ----print=T-----------------------------------------------------------------------
data(carcass)
cm1 <- cmod(~.^., carcass)
cm2 <- cmod(~.^1, data=carcass)

## ----------------------------------------------------------------------------------
testdelete(cm1, ~Meat11:Fat11)
testdelete(cm1, ~Meat12:Fat13)

## ----------------------------------------------------------------------------------
testadd(cm2, ~Meat11:Fat11)
testadd(cm2, ~Meat12:Fat13)

## ----fig=T-------------------------------------------------------------------------
data(carcass)
cm1 <- cmod(~LeanMeat:Meat12:Fat12+LeanMeat:Fat11:Fat12+Fat11:Fat12:Fat13, data=carcass)
plot(cm1)

## ----------------------------------------------------------------------------------
getInEdges(cm1)

## ----------------------------------------------------------------------------------
getInEdges(cm1, type="decomposable")

## ----------------------------------------------------------------------------------
getOutEdges(cm1)

## ----------------------------------------------------------------------------------
getOutEdges(cm1, type="decomposable")

## ----fig=T-------------------------------------------------------------------------
data(carcass)
cm1 <- cmod(~LeanMeat:Meat12:Fat12+LeanMeat:Fat11:Fat12+Fat11:Fat12:Fat13+Fat12:Meat11:Meat13, data=carcass[1:20,])
plot(cm1)

## ----eval=T------------------------------------------------------------------------
in.ed <- getInEdges(cm1)
z <- testInEdges(cm1, edgeList=in.ed)
z

## ----eval=T------------------------------------------------------------------------
z <-testInEdges(cm1, edgeList=in.ed, headlong=T)
z

## ----eval=T------------------------------------------------------------------------
out.ed <- getOutEdges(cm1)
z <- testOutEdges(cm1, edgeList=out.ed)

## ----eval=T------------------------------------------------------------------------
z <- testOutEdges(cm1, edgeList=out.ed, headlong=T)
z

## ----------------------------------------------------------------------------------
dm3 <- dmod(list(c("smoke", "systol"), c("smoke", "mental", "phys")),
            data=reinis)
names(dm3)

## ----------------------------------------------------------------------------------
str(terms(dm3))

## ----------------------------------------------------------------------------------
str(dm3$glistNUM)

## ----------------------------------------------------------------------------------
dm3$varNames

## ----------------------------------------------------------------------------------
str(dm3[c("varNames","conNames","conLevels")])

## ----------------------------------------------------------------------------------
summary(dm1) ## FIXME

## ----------------------------------------------------------------------------------
str(fitted(dm1))
str(dm1$data)

## ----pearson-1,fig=T,include=F, eval=F---------------------------------------------
#  X2 <- (fitted(dm1)-dm1$datainfo$data)/sqrt(fitted(dm1))
#  qqnorm(as.numeric(X2))

