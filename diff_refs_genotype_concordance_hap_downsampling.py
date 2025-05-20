import argparse
import vcfpy

def parse_regions(regions_input):
    """
    Parse user input into a list of regions with chromosomes and positions.
    Expected format: 'chr1:20000-20100' or 'chr1:20000-20100; chr2:30000-31000; chr3:40000-41000'
    """
    regions = []
    for region_str in regions_input.split(";"):
        region_str = region_str.strip()
        if not region_str:
            continue
        try:
            chrom, pos_range = region_str.split(":")
            start_str, end_str = pos_range.split("-")
            regions.append({
                "chrom": chrom.strip(),
                "start": int(start_str),
                "end": int(end_str)
            })
        except Exception as e:
            print(f"Error parsing region '{region_str}': {e}")
    return regions

def process_region(region, truth_vcf_path, call_vcf_path):
    print(f"\n=== Processing region: {region['chrom']}:{region['start']}-{region['end']} ===")
    
    # Parse the truth VCF for the region
    vcf_reader1 = vcfpy.Reader.from_path(truth_vcf_path)
    truth_records = {}
    print("Parsing truth VCF records in region:")
    for record in vcf_reader1:
        if record.CHROM != region["chrom"] or record.POS < region["start"] or record.POS > region["end"]:
            continue
        truth_records[record.POS] = record  # using POS as key
        print(f"Truth record at {record.CHROM}:{record.POS}")
        print("  Reference allele:", record.REF)
        print("  Alternative alleles:", [alt.value for alt in record.ALT])
        for call in record.calls:
            sample_name = call.sample
            genotype = call.data.get("GT")
            print(f"    Sample {sample_name} original genotype: {genotype}")
 
    # Parse the call VCF for the region
    vcf_reader2 = vcfpy.Reader.from_path(call_vcf_path)
    call_records = {}
    print("\nParsing call VCF records in region:")
    for record in vcf_reader2:
        if record.CHROM != region["chrom"] or record.POS < region["start"] or record.POS > region["end"]:
            continue
        call_records[record.POS] = record 
        print(f"Call record at {record.CHROM}:{record.POS}")
        print("  Reference allele:", record.REF)
        print("  Alternative alleles:", [alt.value for alt in record.ALT])
        for call in record.calls:
            sample_name = call.sample
            genotype = call.data.get("GT")
            print(f"    Sample {sample_name} original genotype: {genotype}")
 
    # Compare records and convert genotypes
    print("\nComparing records and converting genotypes:")
    for pos, call_record in call_records.items():
        if pos not in truth_records:
            print(f"Position {pos} found in call VCF but not in truth VCF.")
            continue
 
        truth_record = truth_records[pos]
        # If the reference alleles are identical, use sequence mapping.
        if truth_record.REF == call_record.REF:
            print(f"Position {pos}: Matching reference allele ({truth_record.REF})")
            truth_allele_map = {"0": truth_record.REF}
            for i, alt in enumerate(truth_record.ALT, start=1):
                truth_allele_map[str(i)] = alt.value
 
            call_allele_map = {"0": call_record.REF}
            for i, alt in enumerate(call_record.ALT, start=1):
                call_allele_map[str(i)] = alt.value
            # Define call_allele_map_normal here so it exists in conversion loop.
            call_allele_map_normal = call_allele_map
            mapping_type = "sequence"
        else:
            print(f"Different reference alleles at {pos}:")
            print(f"  Truth REF: {truth_record.REF}")
            print(f"  Call REF:  {call_record.REF}")
            print("Comparing based on differences between alleles and their respective reference alleles.\n")
            
            # Build truth allele mapping: difference = len(allele) - len(truth_REF)
            truth_allele_map = {"0": 0}
            for i, alt in enumerate(truth_record.ALT, start=1):
                truth_allele_map[str(i)] = len(alt.value) - len(truth_record.REF)
            
            # Build call allele mapping using call's own reference:
            call_allele_map_normal = {"0": 0}
            for i, alt in enumerate(call_record.ALT, start=1):
                call_allele_map_normal[str(i)] = len(alt.value) - len(call_record.REF)
            # Also build an adjusted version for the call reference allele:
            call_allele_map_adjusted = dict(call_allele_map_normal)
            call_allele_map_adjusted["0"] = len(call_record.REF) - len(truth_record.REF)
 
            print("  Truth allele difference mapping:", truth_allele_map)
            print("  Call allele difference mapping (normal):", call_allele_map_normal)
            print("  Call allele difference mapping (adjusted):", call_allele_map_adjusted)
            
            mapping_type = "length"
 
        print("  Converting genotypes for this record:")
        for call in call_record.calls:
            sample_name = call.sample
            orig_genotype = call.data.get("GT")
            print(f"    Sample {sample_name} original genotype: {orig_genotype}")
 
            if orig_genotype in ("./.", ".|."):
                print(f"    *Converted genotype remains missing*: {orig_genotype}")
                continue
 
            # Locate the matching truth genotype for possible decision making.
            truth_sample = None
            for tcall in truth_record.calls:
                if tcall.sample == sample_name:
                    truth_sample = tcall
                    break
            truth_genotype = truth_sample.data.get("GT") if truth_sample else ""
 
            delimiter = "/" if "/" in orig_genotype else "|"
            allele_codes = orig_genotype.split(delimiter)
            converted_alleles = []
 
            for allele in allele_codes:
                if allele == ".":
                    converted_alleles.append(".")
                else:
                    match_found = False
                    if mapping_type == "length":
                        # For allele "0", choose between the normal and adjusted mapping based on the truth genotype.
                        if allele == "0":
                            if "0" in truth_genotype:
                                mapping_dict = call_allele_map_normal
                            else:
                                mapping_dict = call_allele_map_adjusted
                        else:
                            mapping_dict = call_allele_map_normal
                        call_allele_difference = mapping_dict.get(allele)
                        for key, truth_difference in truth_allele_map.items():
                            if truth_difference == call_allele_difference:
                                converted_alleles.append(key)
                                match_found = True
                                break
                        if not match_found:
                            print(f"      ALT allele with difference '{call_allele_difference}' from sample {sample_name} not found in the truth VCF at position {pos}.")
                            converted_alleles.append("N/A")
                    else:  # mapping_type == "sequence"
                        mapping_dict = call_allele_map_normal
                        call_allele_value = mapping_dict.get(allele)
                        for key, allele_seq in truth_allele_map.items():
                            if allele_seq == call_allele_value:
                                converted_alleles.append(key)
                                match_found = True
                                break
                        if not match_found:
                            print(f"      ALT allele '{call_allele_value}' from sample {sample_name} not found in the truth VCF at position {pos}.")
                            converted_alleles.append("N/A")
 
            converted_genotype = delimiter.join(converted_alleles)
            print(f"    Sample {sample_name} converted genotype: {converted_genotype}\n")
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="VCF genotype comparison for specified regions")
    parser.add_argument("--truth", type=str, required=True, help="Path to the truth VCF file")
    parser.add_argument("--call", type=str, required=True, help="Path to the call VCF file")
    parser.add_argument("--regions", type=str, required=True, help="Regions to process. Format: 'chr1:20000-20100; chr2:30000-31000'")
    args = parser.parse_args()
 
    regions = parse_regions(args.regions)
    for region in regions:
        process_region(region, args.truth, args.call)
