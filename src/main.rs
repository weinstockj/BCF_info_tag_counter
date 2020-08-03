use std::env;
use std::str;
use rust_htslib::bcf;
use rust_htslib::bcf::{Read};
use log::{info};
use std::collections::HashMap;
use std::fs::File;
use chrono::prelude::*;

fn main() {
    env_logger::init();

    let args: Vec<String>  = env::args().collect();
    let bcf_path = &args[1];

    info!("now running for bcf {}", bcf_path);

    let mut bcf_reader = bcf::Reader::from_path(&bcf_path).unwrap();
    let sample_bytes = bcf_reader.header().samples();
    let samples = sample_bytes.iter().map(|x| String::from_utf8(x.to_vec()).unwrap()).collect::<Vec<String>>(); // convert bytes to strings
    let sample_count = bcf_reader.header().sample_count() as usize;

    info!("sample count: {:?}", sample_count); // print out sample names
    /* info!("samples: {:?}", samples); // print out sample names */

    /* let threshold = 400; */
    let mut count = 0;

    let mut record = bcf_reader.empty_record();

    let mut sample_cosmic_count = HashMap::new();
    let mut variant_count = HashMap::new();

    let mut variant_key;

    loop {
        match bcf_reader.read(&mut record) {
            Ok(true) => (),
            Ok(false) => break,
            Err(_e) => break
        };

    /* for r in bcf_reader.records() { */
        /* if count > threshold { */
        /*     break; */
        /* } */
        /* let mut record = r.unwrap(); */
        let rid = record.rid().unwrap();
        let chrom = bcf_reader.header().rid2name(rid).unwrap().to_owned();
        let chrom = str::from_utf8(&chrom).unwrap().to_owned();
        let pos = record.pos().to_owned() + 1; // rust-htslib returns 0-based values
        let cosmic = record.info(b"COSMIC").flag().unwrap().to_owned();
        let alleles = record.alleles().into_iter().map(|x| str::from_utf8(x).unwrap()).collect::<Vec<&str>>();

        variant_key = format!("{}-{}-{}-{}", chrom, pos, alleles[0], alleles[1]);

        /* info!("variant key: {:?} ", variant_key); */
        
        if cosmic {

            let genotypes = record.genotypes().unwrap();
            let mut alt_across_samples = 0;
            for s in 0..sample_count {
                let mut alt = 0;
                if genotypes.get(s)[0] == bcf::record::GenotypeAllele::Unphased(1) {
                    alt += 1;
                    alt_across_samples += 1;
                }
                if genotypes.get(s)[1] == bcf::record::GenotypeAllele::Unphased(1) {
                    alt += 1;
                    alt_across_samples += 1;
                }

                if sample_cosmic_count.contains_key(&samples[s]) {
                    *sample_cosmic_count.entry(&samples[s]).or_insert(0) += alt;
                } else {
                    sample_cosmic_count.insert(&samples[s], alt);
                }
                /* info!("genotype second haplotype {:?}", genotypes.get(s)[1]); */
            }

            if variant_count.contains_key(&variant_key) {
                *variant_count.entry(variant_key).or_insert(0) += alt_across_samples;
            } else {
                variant_count.insert(variant_key, alt_across_samples);
            }

        }
        /* info!("cosmic: {:?}, ", genotypes); // print out sample names */
        /* info!("cosmic: {:?}, ", cosmic); // print out sample names */
        count = count + 1;
    }

    info!("completed parsing of {} variants", count);
    info!("now writing to disk");

    let local: DateTime<Local> = Local::now();
    let date = local.format("%Y_%m_%d").to_string();

    let mut sample_count_fname = File::create(format!("sample_count_{}.json", date)).unwrap();
    let mut variant_count_fname = File::create(format!("variant_count_{}.json", date)).unwrap();

    serde_json::to_writer(&mut sample_count_fname, &sample_cosmic_count).unwrap();
    serde_json::to_writer(&mut variant_count_fname, &variant_count).unwrap();

    info!("done.");
}
