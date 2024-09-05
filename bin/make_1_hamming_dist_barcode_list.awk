BEGIN {
    FS = ",";  # Set input field separator to comma
    OFS = "\t"; # Set output field separator to tab
    split("NACGT", replacements, ""); # Define the replacement characters without spaces
}

{
    original = $2; # Read the original barcode from the second column
    modified = ""; # Initialize the modified string

    # Iterate over each character of the original barcode
    for (i = 1; i <= length(original); i++) {
        # Extract the current character
        char = substr(original, i, 1);

        # Iterate over the replacement characters
        for (j in replacements) {
            replacement = replacements[j];
            # Skip if the replacement character is the same as the original character
            if (char == replacement) continue;

            # Generate the modified barcode with the current character replaced
            if (modified != "") modified = modified ",";
            modified = modified substr(original, 1, i - 1) replacement substr(original, i + 1);
        }
    }

    # Print the original barcode and the comma-separated list of modified barcodes
    print original, modified;
}