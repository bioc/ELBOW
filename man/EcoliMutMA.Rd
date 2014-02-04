\name{EcoliMutMA}
\docType{data}
\alias{EcoliMutMA}
\title{Escherichia coli str. K-12 substr. MG1655 - EcoliMutMA (GSM576640-45 mutant -- G3.2 mutant of E.coli EcNR1, no isobutanol)}
\description{
	Evolution combined with genomic study elucidates genetic
	bases of isobutanol tolerance in Escherichia coli

	\bold{Platform organism:} Escherichia coli str. K-12 substr. MG1655 \cr
	\bold{Sample organism:}   Escherichia coli \cr
	\bold{Experiment type:}   Expression profiling by array \cr
	\bold{Summary} \cr

	\emph{BACKGROUND:} Isobutanol is a promising next generation biofuel
	with demonstrated high yield microbial production, but the toxicity
	of this molecule reduces fermentation volumetric productivity and final
	titers. Organic solvent tolerance is a complex, multigenic phenotype
	that has been recalcitrant to rational engineering approaches. We apply
	experimental evolution followed by genome resequencing and a gene
	expression study to elucidate genetic bases on adaptation to exogenous
	isobutanol stress. \cr \cr

	\emph{RESULTS:} The adaptations acquired in our evolved lineages exhibit
	antagonistic pleiotropy between minimal and rich medium, and appear to be
	specific to the effects of longer chain alcohols. By examining genotypic
	adaptation in multiple independent lineages, we find evidence of parallel
	evolution in hfq, mdh, acrAB, gatYZABCD, and rph genes. Many isobutanol
	tolerant lineages show reduced rpoS activity, perhaps related to mutations
	in hfq or acrAB. Consistent with the complex, multigenic nature of solvent
	tolerance, we observe adaptations in a diversity of cellular processes.
	Many adaptations appear to involve epistasis between different mutations,
	implying a rugged fitness landscape for isobutanol tolerance. We observe a
	trend of evolution targeting post-transcriptional regulation and high
	centrality nodes of biochemical networks. Collectively, the genotypic
	adaptations we observe suggest mechanisms of adaptation to isobutanol stress
	based on remodelling the cell envelope and surprisingly, stress response
	attenuation. \cr \cr

	\emph{CONCLUSIONS:} We have discovered a set of genotypic adaptations that
	confer increased tolerance to exogenous isobutanol stress. Our results are
	immediately useful to efforts to engineer more isobutanol tolerant host
	strains of E. coli for isobutanol production. We suggest that rpoS and
	post-transcriptional regulators, such as hfq, RNA helicases, and sRNAs
	may be interesting mutagenesis targets for futurue global phenotype engineering. \cr \cr
 	
	\bold{Overall design} \cr
	Two strains (WT strain and G3.2 mutant strain), each with two culture conditions
	(with and without isobutanol in medium). Three biological replicates for each
	strain/culture condition. Twelve samples in total.

	\itemize{
		\item{GSM576634 -- WT E.coli EcNR1, no isobutanol}
		\item{GSM576635 -- WT E.coli EcNR1, no isobutanol}
		\item{GSM576636 -- WT E.coli EcNR1, no isobutanol}
		\item{GSM576637 -- WT E.coli EcNR1, 0.5 \% isobutanol}
		\item{GSM576638 -- WT E.coli EcNR1, 0.5 \% isobutanol}
		\item{GSM576639 -- WT E.coli EcNR1, 0.5 \% isobutanol}
		\item{GSM576640 -- G3.2 mutant of E.coli EcNR1, no isobutanol -- EcoliMutMA}
		\item{GSM576641 -- G3.2 mutant of E.coli EcNR1, no isobutanol}
		\item{GSM576642 -- G3.2 mutant of E.coli EcNR1, no isobutanol}
		\item{GSM576643 -- G3.2 mutant of E.coli EcNR1, 0.5 \% isobutanol}
		\item{GSM576644 -- G3.2 mutant of E.coli EcNR1, 0.5 \% isobutanol}
		\item{GSM576645 -- G3.2 mutant of E.coli EcNR1, 0.5 \% isobutanol}
	}

	\bold{Correponding dataset} \cr
	EcoliMutMA -- GSM576640-45 mutant
}
\usage{EcoliMutMA}
\format{
	A 7 column table defined as such:
	\itemize{
		\item{the first column containing the names of the probes}
		\item{columns 2-4 contain the gene expression values of
			three replicates for the initial conditions of
			the experiment}
		\item{columns 5-7 contain the gene expression values of
			three replicates for the final conditions of
			the experiment}
	}
}
\source{Geo Accession #Series GSE23526}

\references{Minty, J. J., Lesnefsky, A. A., Lin, F., Chen, Y., Zaroff, T. A., Veloso, A. B., Xie, B., McConnell, C. A.,Ward, R. J., Schwartz, D. R., Rouillard, J. M., Gao, Y., Gulari, E., Lin, X.N. (March 2011) \dQuote{Evolution combined with genomic study elucidates genetic bases of isobutanol tolerance in Escherichia coli.} \emph{Microb Cell Fact.} Vol. \bold{10}, p. 18. doi: 10.1186/1475-2859-10-18.}

\keyword{datasets}
