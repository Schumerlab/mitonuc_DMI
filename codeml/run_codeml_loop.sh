#!/bin/bash
## Runs codeml for all mtDNA ND genes, plus NDUFA13 and NDUFS5
## Requires codeml and R to be in $PATH

echo -e "gene\tbir_dNdS\tother_dNdS\tconstant_dNdS\tP_variedrate\tP_nonneutral" > codeml_summary.txt

for gene in $(cat codeml_gene.list); do

    codeml codeml_${gene}_no_branch_rate.ctl
    codeml codeml_${gene}_branch_rate.ctl
    codeml codeml_${gene}_fixedw_branch_rate.ctl

    bir_dNdS=$(grep '(dN/dS) for branches' codeml_${gene}_allowratechange.out | awk '{print $6}')
    other_dNdS=$(grep '(dN/dS) for branches' codeml_${gene}_allowratechange.out | awk '{print $5}')
    constant_dNdS=$(grep 'omega (dN/dS)' codeml_${gene}_samerate.out | awk '{print $4}')
    nochange=$(grep 'lnL' codeml_${gene}_samerate.out | awk '{print $5}')
    allowchange=$(grep 'lnL' codeml_${gene}_allowratechange.out | awk '{print $5}')
    fixedw=$(grep 'lnL' codeml_${gene}_fixedw_allowratechange.out | awk '{print $5}')
    
    P_variedrate=$(R -e "pchisq(2*($allowchange - $nochange), df=1, lower.tail=FALSE)" | grep '\[1\]' | awk '{print $2}')
    P_nonneutral=$(R -e "pchisq(2*($allowchange - $fixedw), df=1, lower.tail=FALSE)" | grep '\[1\]' | awk '{print $2}')

    echo -e "${gene}\t${bir_dNdS}\t${other_dNdS}\t${constant_dNdS}\t${P_variedrate}\t${P_nonneutral}" >> codeml_summary.txt

    done
