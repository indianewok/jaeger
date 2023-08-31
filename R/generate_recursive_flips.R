Rcpp::cppFunction('
  #include <Rcpp.h>
  #include <bitset>
  #include <set>
  using namespace Rcpp;

  // Recursive function to flip bits
  void recursive_flip_bits(int64_t original_value, int64_t value, int depth, int max_depth, int num_bits, std::set<int64_t>& results) {
    if (depth >= max_depth) return; // Base case

    for (int index = 0; index < num_bits; ++index) { // Loop through the specified number of bits
      int64_t mask = 1LL << (num_bits - 1 - index); // Shift the 1 to the target position
      int64_t result = value ^ mask; // Use XOR to flip the specified bit
      if (result != original_value) { // Exclude the original value
        results.insert(result); // Add the result to the set
        recursive_flip_bits(original_value, result, depth + 1, max_depth, num_bits, results); // Recursive call
      }
    }
  }

  // [[Rcpp::export]]
  NumericVector generate_recursive_flips(SEXP value, int max_depth, int barcode_length) {
    NumericVector numeric_value(value); // Convert SEXP to NumericVector
    int64_t int64_value = static_cast<int64_t>(numeric_value[0]); // Cast the first element to int64_t
    int num_bits = 2 * barcode_length; // Number of bits corresponding to the barcode length

    std::set<int64_t> unique_results; // Set to store unique flipped results
    recursive_flip_bits(int64_value, int64_value, 0, max_depth, num_bits, unique_results); // Call the recursive function

    return wrap(std::vector<int64_t>(unique_results.begin(), unique_results.end())); // Wrap the set as integer64
}')
