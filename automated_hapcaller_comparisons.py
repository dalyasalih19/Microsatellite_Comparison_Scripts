import sys
import argparse
import vcfpy

def parse_regions(regions_input):
  """Parse user input into a list of region with chromosomes and positions. Expected format: 'chr1:20000-20100' or 'chr1:20000-20100; chr2:30000-31000; chr3:40000-41000'"""
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

def main():
  parser = argparse.ArgumentParser(description="Automated Hapcaller Comparisons")
  parser.add_argument("-v", "--vcf", required=True, help="Path to the VCF file")
  parser.add_argument("-r", "--regions", required=True, help="Region(s) in the format: 'chr1:20000-20100' or 'chr1:20000-20100; chr2:30000-31000; chr3:40000-41000'")
  args = parser.parse_args()

  vcf_file = args.vcf.strip()
  if not vcf_file:
    sys.exit("Error: VCF file path is missing.")

  regions_input = args.regions.strip()
  if not regions_input:
    sys.exit("Error: Regions input is missing.")

  regions = parse_regions(regions_input)
  if not regions:
    sys.exit("Error: No valid regions provided.")

  # Open VCF
  try:
    vcf_reader = vcfpy.Reader.from_path(vcf_file)
  except Exception as e:
    sys.exit(f"Error opening VCF file: {e}")
    
  print("VCF header:")
  print(str(vcf_reader.header))
    
  print("VCF contig lines:")
  for contig_line in vcf_reader.header.get_lines("contig"):
    print(contig_line)

  # Get available chromosomes from the VCF header.
  available_chroms = [line.id for line in vcf_reader.header.get_lines("contig")]
  for r in available_chroms:
    print("Avaliable chromosomes: ", available_chroms)
    print("I am within this loop")
  print("Just before invalid regions")
  invalid_regions = [r for r in regions if r["chrom"] not in available_chroms]
  if invalid_regions:
    for r in invalid_regions:
      print(f"Warning: Chromosome {r['chrom']} not found in VCF file. Skipping region {r['chrom']}:{r['start']}-{r['end']}")
  # Filter out invalid regions.
  regions = [r for r in regions if r["chrom"] in available_chroms]
  if not regions:
    sys.exit("Error: None of the specified regions are present in the VCF.")

  # Read all records into memory (assumes VCF is of manageable size)
  all_records = list(vcf_reader)

  # Process each valid region.
  for region in regions:
    region_records = []
    region_chrom = region["chrom"]
    region_start = region["start"]
    region_end = region["end"]
    print(f"\nProcessing region {region_chrom}:{region_start}-{region_end}\n")
    
    # Filter records for this region.
    region_records = [record for record in all_records if record.CHROM == region_chrom and region_start <= record.POS <= region_end]
    
    if not region_records:
      print(f"Error: No records found in region {region_chrom}:{region_start}-{region_end}")
      continue

        # Count the number of indels in the region.
    indel_count = 0
    for record in region_records:
      # Check if any alt allele differs in length from the reference.
      for alt in record.ALT:
        if alt is not None and len(record.REF) != len(str(alt.value)):
          indel_count += 1
          break  # Only count this record once.
    if indel_count > 1:
      print(f"\n***COMMENT: Complex region with {indel_count} indels in the VCF***\n")

    for record in region_records:
      print("Processing record:", record.CHROM, record.POS)
      
      # Get reference and alternate alleles.
      ref = record.REF
      alts = [alt.value for alt in record.ALT] if record.ALT else []
      print("Reference allele:", ref, "Alternative alleles:", alts)
      
      # Build the temp array: first element is length(ref), then lengths of each alt.
      temp_array = [len(ref)]
      for alt in alts:
        if alt is not None:
          temp_array.append(len(str(alt)))
      print("Temp array of allele lengths:", temp_array)
      
      # Find the smallest element in the temp array.
      smallest = min(temp_array)
      print("Smallest element:", smallest)
      
      # Create a normalized (subtracted) array.
      sub_array = [x - smallest for x in temp_array]
      print("Subtracted array:", sub_array)
      
      # Process genotype calls for each sample in this record.
      print("The Original Genotypes:")
      for call in record.calls:
        sample_name = call.sample
        genotype = call.data.get("GT")
        print("Processing sample:", sample_name, "genotype:", genotype)
        # Replace phased delimiter "|" with "/" and split.
        genotype_clean = genotype.replace("|", "/")
        alleles = genotype_clean.split("/")
        genotype_lengths = []
        for allele in alleles:
          try:
            allele_index = int(allele)
            if allele_index < len(sub_array):
              genotype_lengths.append(sub_array[allele_index])
            else:
              genotype_lengths.append("NA")
          except ValueError:
            genotype_lengths.append("NA")
        # Join the new genotype lengths with "/" as delimiter.
        new_genotype = "/".join(str(x) for x in genotype_lengths)
        print("Sample:", sample_name, "| Original genotype:", genotype, "| New genotype lengths:", genotype_lengths)

if __name__ == "__main__":
    main()
