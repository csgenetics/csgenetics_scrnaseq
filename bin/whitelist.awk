BEGIN {
    OFS="\t"; # Set output field separator to tab
}

{
    original = $1; # Read the original barcode from the first column
    modified = ""; # Initialize the modified string

    # Iterate over each character of the original barcode
    for (i = 1; i <= length(original); i++) {
        # Extract the current character
        char = substr(original, i, 1);

        # Iterate over the replacement characters
        for (replacement in replacements) {
            # Skip if the replacement character is the same as the original character
            if (char == replacements[replacement]) continue;

            # Generate the modified barcode with the current character replaced
            if (modified != "") modified = modified ",";
            modified = modified substr(original, 1, i - 1) replacements[replacement] substr(original, i + 1);
        }
    }

    # Print the original barcode and the comma-separated list of modified barcodes
    print original, modified;
}

# Define the replacement characters
BEGIN {
    split("N A C T G", replacements);
}