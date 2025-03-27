

## select LoF and damaging missense varaints for further gene-based analysis
perl selectVariants.pl

## transver the variant annotation file into a gene-case matrix for association and permutation
perl generateMatrix.pl

## perform association
for pheno in (ALL EA MIX FF MA); do
	geneBurden.R $pheno
done

## perform permutation
thread=10
for pheno in (ALL EA MIX FF MA); do
	perl runPerm1k.pl $pheno $thread
done

## plot distributions of P values
for pheno in (ALL EA MIX FF MA); do
	QQplot.R $pheno 
done