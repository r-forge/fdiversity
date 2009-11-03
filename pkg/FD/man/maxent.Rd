\encoding{latin1}
\name{maxent}
\Rdversion{1.1}
\alias{maxent}
\concept{maximum entropy}
\concept{statistical mechanics}
\title{
Estimating Probabilities via Maximum Entropy: Improved Iterative Scaling
}
\description{
\code{maxent} returns the probabilities that maximize the entropy conditional on a series of constraints that are linear in the features. It relies on the Improved Iterative Scaling algorithm of Della Pietra et al. (1997). It has been used to obtain the predicted relative abundances of a set of species given the trait values of each species and the community-aggregated trait values at a site (Shipley et al. 2006; Shipley 2009; Sonnier et al. 2009).
}

\usage{
maxent(c.means, c.mat, prior, tol = 1e-08, lambda = F)
}

\arguments{
  \item{c.means}{
vector of the constraint means (e.g. community-aggregated trait values at each site)
}
  \item{c.mat}{
matrix of constraints (rows) for each state (columns), e.g. a traits (rows) by species (columns) matrix
}
  \item{prior}{
vector of prior probabilities. Can be missing, in which case a uniform distribution is assumed (i.e. maximally uninformative prior).
}
  \item{tol}{
tolerance threshold to determine convergence. See \sQuote{details} section. 
}
  \item{lambda}{
Logical. Should \eqn{\lambda}{lambda}-values be returned? 
}
}
\details{
The biological model of community assembly through trait-based habitat
filtering (Keddy 1992) has been translated mathematically
via a maximum entropy (maxent) model by Shipley et al. (2006) and
Shipley (2009). A maxent model contains three components: (i) a set
of possible states, (ii) a set of macroscopic empirical constraints,
and (iii) a prior probability distribution \eqn{\mathbf{q}=[q_j]}{q = [qj]}.

In the context of community assembly, states are species, macroscopic
empirical constraints are community-aggregated traits, and prior probabilities
\eqn{\mathbf{q}}{q} are the relative abundances of species of the regional
pool (Shipley et al. 2006, Shipley 2009). By default, these prior
probabilities \eqn{\mathbf{q}}{q} are uniformally distributed (a maximally uninformative
prior), but can be specificied otherwise (Shipley 2009, Sonnier et
al. 2009).

To facilitate the link between the biological model and the mathematical
model, in the following description of the algorithm states are species and constraints are traits.

If \eqn{\mathbf{q}}{q} is not specified, set \eqn{p_{j}=1/S}{pj = 1 / S} for each of the
\eqn{S}{S} species (i.e. a uniform distribution), where \eqn{p_{j}}{pj} is the
probability of species \eqn{j}{j}, otherwise \eqn{p_{j}=q_{j}}{pj = qj}.

Calulate a vector \eqn{\mathbf{c=\left[\mathrm{\mathit{c_{i}}}\right]}=\{c_{1},\; c_{2},\;\ldots,\; c_{T}\}}{c = [ci] = {c1, c2, ..., cT}},
where \eqn{c_{i}={\displaystyle \sum_{j=1}^{S}t_{ij}}}{ci = sum(tij)}; i.e. each \eqn{c_{i}}{ci}
is the sum of the values of trait \eqn{i}{i} over all species, and \eqn{T}{T}
is the number of traits.

Repeat for each iteration \eqn{k}{k} until convergence:

1. For each trait \eqn{t_{i}}{ti} (i.e. row of the constraint matrix) calculate:

\deqn{
\gamma_{i}(k)=ln\left(\frac{\bar{t}_{i}}{{\displaystyle \sum_{j=1}^{S}\left(p_{j}(k)\; t_{ij}\right)}}\right)\left(\frac{1}{c_{i}}\right)}{%
gamma_i(k) = ln(t.mean_i / [sum(pj(k) tij)]) (1 / ci)}

This is simply the natural log of the known community-aggregated
trait value to the calculated community-aggregated trait value at
this step in the iteration, given the current values of the probabilities.
The whole thing is divided by the sum of the known values of the trait
over all species.

2. Calculate the normalization term \eqn{Z}{Z}:

\deqn{Z(k)=\left({\displaystyle \sum_{j=1}^{S}p_{j}(k)\; e^{\left({\displaystyle \sum_{i=1}^{T}\gamma_{i}(k)}\; t_{ij}\right)}}\right)}{%
Z(k) = sum(pj(k) e^(gamma_i(k) tij) ) }

3. Calculate the new probabilities \eqn{p_{j}}{pj} of each species at iteration \eqn{k+1}{k+1}:

\deqn{p_{j}(k+1)=\frac{{\displaystyle p_{j}(k)\; e^{\left({\displaystyle \sum_{i=1}^{T}\gamma_{i}(k)}\; t_{ij}\right)}}}{Z(k)}}{%
pj(k+1) = [pj(k) e^(gamma_i(k) tij)/ Z]}

4. If \eqn{|max\left(p\left(k+1\right)-p\left(k\right)\right)|\leq}{|max(pj(k+1) - pj(k))| <=} tolerance threshold (i.e. argument \code{tol}) then stop, else repeat steps 1 to 3.

When convergence is achieved then the resulting probabilities (\eqn{\hat{p}_{j}}{pj_hat})
are those that are as close as possible to \eqn{q_j}{qj} while simultaneously maximize the entropy conditional on
the community-aggregated traits and that maximize the likelihood of
the Gibbs distribution:

\deqn{\hat{p}_{j}=\frac{q_{j}e^{\left({\displaystyle -}{\displaystyle \sum_{i=1}^{T}\lambda_{i}t_{ij}}\right)}}{{\displaystyle \sum_{j=1}^{S}q_{j}}e^{\left({\displaystyle -}{\displaystyle \sum_{i=1}^{T}\lambda_{i}t_{ij}}\right)}}=\frac{q_{j}e^{\left({\displaystyle -}{\displaystyle \sum_{i=1}^{T}\lambda_{i}t_{ij}}\right)}}{Z}}{%
Note: equation not shown in HTML help file: please refer to PDF manual.}

This means that one can solve for the Langrange multipliers (i.e.
weights on the traits, \eqn{\lambda_{i}}{lamda_i}) by solving the linear system
of equations:

\deqn{\left(\begin{array}{c}
ln\left(\hat{p}_{1}\right)\\
ln\left(\hat{p}_{2}\right)\\
\vdots\\
ln\left(\hat{p}_{S}\right)\end{array}\right)=\left(\lambda_{1},\;\lambda_{2},\;\ldots,\;\lambda_{T}\right)\left[\begin{array}{cccc}
t_{11} & t_{12} & \ldots & t_{1S}-ln(Z)\\
t_{21} & t_{22} & \vdots & t_{2S}-ln(Z)\\
\vdots & \vdots & \vdots & \vdots\\
t_{T1} & t_{T2} & \ldots & t_{TS}-ln(Z)\end{array}\right]-ln(Z)}{%
Note: equation not shown in HTML help file: please refer to PDF manual.}

This system of linear equations has \eqn{T+1}{T+1} unknowns (the \eqn{T}{T} values
of \eqn{\lambda}{lambda} plus \eqn{ln(Z)}{ln(Z)}) and \eqn{S}{S} equations. So long as the number
of traits is less than \eqn{S-1}{S-1}, this system is soluble. In fact, the
solution is the well-known least squares regression: simply regress
the values \eqn{ln(\hat{p}_{j})}{ln(pj_hat} of each species on the trait values
of each species in a multiple regression.

The intercept is the value of \eqn{ln(Z)}{ln(Z)} and the slopes are the values
of \eqn{\lambda_{i}}{lambda_i} and these slopes (Lagrange multipliers) measure
by how much the \eqn{ln(\hat{p}_{j})}{ln(pj_hat)}, i.e. the \eqn{ln}{ln}(relative abundances),
changes as the value of the trait changes.
}

\value{
  \item{prob}{vector of predicted probabilities}
  \item{moments}{vector of final moments}
  \item{entropy}{Shannon entropy of \code{prob}}
  \item{iter}{number of iterations required to reach convergence}
  \item{lambda}{\eqn{\lambda}{lambda}-values, only returned if \code{lambda = T}}
}
\references{
Della Pietra, S., V. Della Pietra, and J. Lafferty (1997) Inducing features of random fields. \emph{IEEE Transactions Pattern Analysis and Machine Intelligence} \bold{19}:1-13.

Keddy, P. A. (1992) Assembly and response rules: two goals for predictive community ecology. \emph{Journal of Vegetation Science} \bold{3}:157-164.

Shipley, B., D. Vile, and \enc{�}{E}. Garnier (2006) From plant traits to plant communities: a statistical mechanistic approach to biodiversity. \emph{Science} \bold{314}: 812--814.

Shipley, B. (2009) From Plant Traits to Vegetation Structure: Chance and Selection in the Assembly of Ecological Communities. Cambridge University Press, Cambridge, UK. 290 pages.

Sonnier, G., Shipley, B., and M. L. Navas. 2009. Plant traits, species pools and the prediction of relative abundance in plant communities: a maximum entropy approach. \emph{Journal of Vegetation Science} \bold{in press}.
}
\author{
Bill Shipley \email{bill.shipley@usherbrooke.ca}

\url{http://pages.usherbrooke.ca/jshipley/recherche/}

Ported to \pkg{FD} by Etienne \enc{Lalibert�}{Laliberte}.
}

\seealso{
\code{\link{functcomp}} to compute community-aggregated traits.
}
\examples{
# an unbiased 6-sided dice, with mean = 3.5
# what is the probability associated with each side,
# given this constraint?
maxent(3.5, 1:6)

# a biased 6-sided dice, with mean = 4
maxent(4, 1:6)


# example with tussock dataset
traits <- tussock$trait[, c(2:7, 11)] # use only continuous traits
traits <- na.omit(traits) # remove 2 species with NA's
labs <- rownames(traits)
# standardise traits to zero mean and unit variance
traits <- apply(traits, 2, scale, center = T)
# add abs(min(traits)) to make sure no negative values in traits
traits <- traits + abs(min(traits))
rownames(traits) <- labs
abun <- tussock$abun[, rownames(traits)] # abundance matrix
abun <- t(apply(abun, 1, function(x) x / sum(x) )) # relative abundances
agg <- functcomp(traits, abun) # community-aggregated traits
traits <- t(traits) # transpose matrix

# run maxent on site 1 (first row of abun), all species
pred.abun <- maxent(as.numeric(agg[1,]), traits)

# plot results, note poor fit
y <- pred.abun$prob
x <- as.numeric(abun[1,])
plot(y ~ x, xlim = c(0, 1), ylim = c(0, 1), xlab = "observed abundances",
ylab = "predicted abundances", main = "site 1, all 51 species\nall 7 traits")
abline(0, 1, col = "grey25")
}

\keyword{distribution}
\keyword{math}
\keyword{models}