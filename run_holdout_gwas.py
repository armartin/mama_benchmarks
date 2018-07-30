from __future__ import print_function
import sys
# from gnomad_hail.utils.slack_utils import *
from hail import *
from pprint import pprint
import argparse
hc = HailContext()

def main(args):
    print('Read/assemble phenotypes')
    """
    need to get these from a bunch of split up tables
    if it's not possible to get these from hail easily, potentially read them in with pandas, 
    where you can scan the header before deciding which phenotypes to read
    """
    sex = 'both_sexes'
    kt_phenotype = hc.import_table('gs://phenotype_31063/ukb31063.phesant_phenotypes.both_sexes.tsv.bgz',
                                   key='userId', quote='"', impute=True, types={'userId': TString()}, missing='')
    phenos = {'userId': 's', '50': 'height', '21001': 'bmi', '48': 'waist', '49': 'hip', '23099': 'body_fat', '4080': 'systolic_bp',
              '4079': 'diastolic_bp', '1558': 'alcohol', '1488': 'tea'} # '3456': 'cigarettes',

    kt_phenotype = kt_phenotype.select(phenos.keys())
    kt_phenotype = kt_phenotype.annotate('whr = `48`/`49`')
    kt_phenotype = kt_phenotype.rename(phenos)

    gwas_inds = hc.import_table('gs://prs_ukbb/holdout_gwas/ukb31063.gwas_samples.gwas_vs_holdout.txt', key='s',
                                impute=True, types={'s': TString()})

    gwas_inds = gwas_inds.filter('in_gwas')
    kt_phenotype = kt_phenotype.join(gwas_inds)


    print('Read covariates, QC sites')
    vds_variants = hc.read('gs://phenotype_31063/hail/ukb31063.gwas_variants.autosomes.with_qc_annotations.vds')
    kt_covariates = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_covariates.{}.kt'.format(sex))
    # 655633986626-compute@developer.gserviceaccount.com does not have storage.buckets.get access to phenotype_31063.

    print('Reading UKBB imputed variants')
    import_expr = '{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}'

    covariate_expr = ['sa.covariates.age',
                      'sa.covariates.age_squared',
                      'sa.covariates.isFemale',
                      'sa.covariates.age_isFemale',
                      'sa.covariates.age_squared_isFemale'] + \
                     ['sa.covariates.PC{:}'.format(i) for i in xrange(1, 21)]

    pheno_codes = list(set(kt_phenotype.columns) - set(['s', 'waist', 'hip', 'in_gwas']))


    vds = (
    hc.import_bgen('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{}_v3.bgen'.format(import_expr),
                   sample_file='gs://phenotype_31063/ukb31063.{}.sample'.format('autosomes'),
                   tolerance=0.2)
        .annotate_variants_vds(vds_variants, expr='va.AF = vds.qc.AF, va.info = vds.info')
        .filter_variants_expr('isDefined(va.AF)', keep=True)
        .annotate_samples_table(kt_covariates, root='sa.covariates')
        .annotate_samples_table(kt_phenotype, root='sa.phenotypes')
    )

    vds.variants_table().write('gs://prs_ukbb/holdout_gwas/pheno_31063_holdout_gwas.kt', args.overwrite)
    #kt_results = hc.read_table('gs://prs_ukbb/holdout_gwas/pheno_31063_holdout_gwas.kt')

    print('Run regressions per phenotype')
    for i in range(len(pheno_codes)):
        my_vds = vds.linreg3(ys=['sa.phenotypes.`{}`'.format(pheno_codes[i])],
                             covariates=covariate_expr,
                             use_dosages=True)
        # my_vds = my_vds.annotate_variants_expr('va.results' = ['va.linreg'])
        kt_results = my_vds.variants_table()
        kt_export = kt_results.annotate(['chr = v.contig',
                                         'pos = v.start',
                                         'ref = v.ref',
                                         'alt = v.alt',
                                         'rsid = va.rsid',
                                         'nCompleteSamples = va.linreg.nCompleteSamples',
                                         'AC = va.linreg.AC',
                                         'ytx = va.linreg.ytx[0]',
                                         'beta = va.linreg.beta[0]',
                                         'se = va.linreg.se[0]',
                                         'tstat = va.linreg.tstat[0]',
                                         'pval = va.linreg.pval[0]'])

        kt_export2 = kt_export.drop(['v', 'va'])
        kt_export2.export('gs://prs_ukbb/holdout_gwas/pheno_31063_holdout_gwas_{}.txt.gz'.format(pheno_codes[i]))


if __name__ == '__main__':
    print('Starting run')
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    main(args)
    #  try_slack('@armartin', main, args)

