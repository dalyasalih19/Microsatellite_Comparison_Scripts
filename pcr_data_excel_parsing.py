import csv

# CSV filename and location
filename = 'fill in'

# Read all rows from the CSV (first column is sample name, then each pair of columns represents allele lengths for one region)
samples = []
with open(filename, 'r') as file:
    reader = csv.reader(file)
    header = next(reader)
    for row in reader:
        samples.append(row)

# Determine number of regions (first column is sample names, then every 2 columns represent a region)
num_regions = (len(header) - 1) // 2

# Process each region individually.
for region in range(num_regions):
    # Determine the column indices for this region
    col1 = 2 * region + 1
    col2 = 2 * region + 2
    # Define a region label using the header names
    region_name = f"{header[col1]} & {header[col2]}"
    
    # Collect all allele values for the current region
    all_allele_values = []
    for row in samples:
        try:
            allele1 = float(row[col1])
            allele2 = float(row[col2])
            all_allele_values.extend([allele1, allele2])
        except ValueError:
            continue

    # Create a sorted list of unique allele values.
    temp_array = sorted(set(all_allele_values))
    print(f"\nRegion {region+1} ({region_name}):")
    print("Temp array:")
    for idx, value in enumerate(temp_array):
        print(f"{idx}\t{value}")

    # Find the smallest allele value
    smallest = temp_array[0] if temp_array else None
    print("Smallest element:", smallest)

    # Create the normalized array by subtracting the smallest value from each allele value.
    sub_array = [value - smallest for value in temp_array] if smallest is not None else []
    print("Subtracted array:")
    for idx, value in enumerate(sub_array):
        print(f"{idx}\t{value}")

    # Build a mapping from original allele length to its normalized value.
    norm_mapping = {orig: norm for orig, norm in zip(temp_array, sub_array)}

    # Process each sample and "genotype" them using the normalized mapping
    print("\nGenotypes:")
    for row in samples:
        sample_name = row[0]
        try:
            allele1 = float(row[col1])
            allele2 = float(row[col2])
        except ValueError:
            allele1, allele2 = None, None

        if allele1 is not None and allele2 is not None:
            new_val1 = norm_mapping.get(allele1, "NA")
            new_val2 = norm_mapping.get(allele2, "NA")
            genotype = f"{new_val1}/{new_val2}"
        else:
            genotype = "NA/NA"

        # Print sample name, region, and genotype
        print(f"{sample_name}\t({region_name})\t{genotype}")
    
    # Print a separator line for clarity between regions.
    print("\n" + "="*40)
