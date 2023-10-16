// [[Rcpp::plugins(openmp)]]  // Enable OpenMP in Rcpp
#include "edlib.h"
#include "RcppInt64.h"
#include <Rcpp.h>
#include <omp.h>  // Include OpenMP header for multi-threading
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdint.h>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <regex>
#include <unordered_map>
#include <cstdint>
#include <unordered_set>
#include <numeric>
#include <cmath>

using namespace std;
using namespace Rcpp;

std::string revcomp(const std::string& sequence) {
  std::unordered_map<char, char> complement {
    {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'},
    {'a', 't'}, {'t', 'a'}, {'c', 'g'}, {'g', 'c'},
    {'N', 'N'}, {'n', 'n'}
  };
  std::string reversed(sequence.rbegin(), sequence.rend());
  for (char& nucleotide : reversed) {
    nucleotide = complement[nucleotide];
  }
  return reversed;
}

int dna_hamming_distance(const std::string& a, const std::string& b) {
  int count = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i] != b[i]) {
      ++count;
    }
  }
  return count;
}

std::vector<int64_t> barcodes_to_bits_cpp(const std::vector<std::string>& barcodes) {
  int n = barcodes.size();
  std::vector<int64_t> results(n);
  for(int i = 0; i < n; ++i) {
    const std::string& barcode = barcodes[i];
    int64_t result = 0;
    for (char c : barcode) {
      result <<= 2;
      switch (c) {
      case 'A': break;
      case 'C': result |= 1; break;
      case 'T': result |= 2; break;
      case 'G': result |= 3; break;
      }
    }
    results[i] = result;
  }
  return results;
}

std::vector<std::string> bits_to_barcodes_cpp(const std::vector<int64_t>& input, int barcode_length = 16) {
  int n = input.size();
  std::vector<std::string> results(n);
  for(int i = 0; i < n; ++i) {
    int64_t int64_code = input[i];
    std::string result;
    for (int j = 0; j < barcode_length; ++j) {
      switch (int64_code & 3) {
      case 0: result.insert(result.begin(), 'A'); break;
      case 1: result.insert(result.begin(), 'C'); break;
      case 2: result.insert(result.begin(), 'T'); break;
      case 3: result.insert(result.begin(), 'G'); break;
      }
      int64_code >>= 2;
    }
    results[i] = result;
  }
  return results;
}

void recursive_flip_bits_cpp(int64_t original_value, int64_t value, int depth, int max_iterations,
  int num_bits, const std::unordered_set<int64_t>& whitelist,
  int& min_hamming_distance, std::set<int64_t>& min_hamming_results) {
  if (depth >= max_iterations) return;
  for (int index = 0; index < num_bits; ++index) {
    int64_t mask = 1LL << (num_bits - 1 - index);
    int64_t result = value ^ mask;
    if (result != original_value) {
      if (whitelist.find(result) != whitelist.end()) {
        int hamming_dist = dna_hamming_distance(bits_to_barcodes_cpp({original_value})[0], bits_to_barcodes_cpp({result})[0]);
        if (hamming_dist < min_hamming_distance) {
          min_hamming_distance = hamming_dist;
          min_hamming_results.clear();
          min_hamming_results.insert(result);
        } else if (hamming_dist == min_hamming_distance) {
          min_hamming_results.insert(result);
        }
      }
      recursive_flip_bits_cpp(original_value, result, depth + 1, max_iterations, num_bits, whitelist, min_hamming_distance, min_hamming_results);
    }
  }
}

void recursive_flip_bits_list(int64_t original_value, int64_t value, int depth, int max_mutations,
  int num_bits, std::set<int64_t>& results, std::unordered_set<int64_t>& whitelist) {
  if (depth >= max_mutations) return;
  for (int index = 0; index < num_bits; ++index) {
    int64_t mask = 1LL << (num_bits - 1 - index);
    int64_t result = value ^ mask;
    if (result != original_value) {
      if (whitelist.find(result) != whitelist.end()) {
        results.insert(result);
      }
      recursive_flip_bits_list(original_value, result, depth + 1, max_mutations, num_bits, results, whitelist);
    }
  }
}

std::map<int64_t, std::vector<std::pair<int64_t, int>>> 
  mutate_and_check_list_cpp(const std::vector<int64_t>& barcodes, 
    const std::unordered_set<int64_t>& whitelist, 
    int max_iterations = 4, int max_mutations = 4, 
    int barcode_length = 16) {
    std::map<int64_t, std::vector<std::pair<int64_t, int>>> results;
    int num_bits = 2 * barcode_length;
    for (const auto& int64_barcode : barcodes) {
      int min_hamming_distance = max_mutations; // Initialize with max_mutations
      std::set<int64_t> min_hamming_results;
      recursive_flip_bits_cpp(int64_barcode, int64_barcode, 0, max_iterations, num_bits, 
        whitelist, min_hamming_distance, min_hamming_results);
      if (!min_hamming_results.empty()) {
        std::vector<std::pair<int64_t, int>> mutated_barcodes;
        for (const auto& mutated : min_hamming_results) {
          mutated_barcodes.emplace_back(mutated, min_hamming_distance);
        }
        results[int64_barcode] = mutated_barcodes;
      }
    }
    return results;
  }

std::unordered_map<std::string, std::vector<int>> sig_parse_ext(const std::string& signature) {
  std::unordered_map<std::string, std::vector<int>> parsedSignature;
  std::stringstream ss(signature);
  std::string token;
  while (std::getline(ss, token, '|')) {
    size_t pos = token.find(":");
    std::string id = token.substr(0, pos);
    std::string coords = token.substr(pos + 1);
    std::stringstream coord_ss(coords);
    std::string coord;
    std::vector<int> coord_vec;
    while (std::getline(coord_ss, coord, ':')) {
      coord_vec.push_back(std::stoi(coord));
    }
    parsedSignature[id] = coord_vec;
  }
  return parsedSignature;
}

std::vector<int> find_mismatch_positions(const std::string& original, const std::string& mutated) {
  std::vector<int> mismatch_positions;
  for (size_t i = 0; i < original.size(); ++i) {
    if (original[i] != mutated[i]) {
      mismatch_positions.push_back(i + 1);  // Add 1 to convert from 0-based to 1-based index
    }
  }
  return mismatch_positions;
}

std::pair<std::vector<std::string>, double> best_matches_jaccard_distance_qgram(
    const std::string& original,
  const std::vector<std::string>& candidates, int q) {
  double min_jaccard_distance = 1.0;
  std::vector<std::string> best_matches;
  for (const auto& candidate : candidates) {
    std::unordered_set<std::string> setA, setB;
    for (int i = 0; i <= static_cast<int>(original.size()) - q; ++i) {
      setA.insert(original.substr(i, q));
    }
    for (int i = 0; i <= static_cast<int>(candidate.size()) - q; ++i) {
      setB.insert(candidate.substr(i, q));
    }
    std::unordered_set<std::string> intersection;
    for (const auto& item : setA) {
      if (setB.find(item) != setB.end()) {
        intersection.insert(item);
      }
    }
    double jaccard_index = static_cast<double>(intersection.size()) / (setA.size() + setB.size() - intersection.size());
    double jaccard_distance = 1 - jaccard_index;
    if (jaccard_distance < min_jaccard_distance) {
      min_jaccard_distance = jaccard_distance;
      best_matches.clear();
      best_matches.push_back(candidate);
    } else if (jaccard_distance == min_jaccard_distance) {
      best_matches.push_back(candidate);
    }
  }
  return std::make_pair(best_matches, min_jaccard_distance);
}

// [[Rcpp::export]]
Rcpp::StringVector revcompR(Rcpp::StringVector sequences) {
  int n = sequences.size();
  Rcpp::StringVector results(n);
  int j = 0; // Counter for results
  for (int i = 0; i < n; ++i) {
    if (sequences[i] == NA_STRING) {
      results[j] = NA_STRING;
      continue;
    }
    results[j] = revcomp(Rcpp::as<std::string>(sequences[i]));
    ++j;
  }
  return results;
}

// [[Rcpp::export]]
LogicalVector check_whitelist(DataFrame r_whitelist_df, NumericVector mutations) {
  NumericVector whitelist_bcs = r_whitelist_df["whitelist_bcs"];
  std::unordered_set<int64_t> whitelist_set;
  for (int i = 0; i < whitelist_bcs.size(); ++i) {
    int64_t whitelist_entry = static_cast<int64_t>(whitelist_bcs[i]);
    whitelist_set.insert(whitelist_entry);
  }
  LogicalVector in_whitelist(mutations.size());
  for (int i = 0; i < mutations.size(); ++i) {
    int64_t mutated_barcode = static_cast<int64_t>(mutations[i]);
    in_whitelist[i] = (whitelist_set.find(mutated_barcode) != whitelist_set.end());
  }
  return in_whitelist;
}

// [[Rcpp::export]]
Rcpp::NumericVector barcodes_to_bits(Rcpp::StringVector barcodes) {
  int n = barcodes.size();
  std::vector<int64_t> results(n);  // Changed to int64_t
  for(int i = 0; i < n; ++i) {
    std::string barcode = Rcpp::as<std::string>(barcodes[i]);
    int64_t result = 0;  // Changed to int64_t
    for (char &c : barcode) {
      result <<= 2;
      switch (c) {
      case 'A': break;
      case 'C': result |= 1; break; 
      case 'T': result |= 2; break;
      case 'G': result |= 3; break;
      }
    } 
    results[i] = result;
  } 
  return Rcpp::toInteger64(results);
}

// [[Rcpp::export]]
Rcpp::StringVector bits_to_barcodes(Rcpp::NumericVector input, int barcode_length = 16, bool verbose = false) {
  // Convert the NumericVector to std::vector<int64_t>
  std::vector<int64_t> cpp_bits = Rcpp::fromInteger64(input);
  int n = cpp_bits.size();
  Rcpp::StringVector results(n);
  for(int i = 0; i < n; ++i) {
    int64_t int64_code = cpp_bits[i];  // Changed to int64_t and used cpp_bits
    std::string result;
    for (int j = 0; j < barcode_length; ++j) {
      if (verbose) {
        Rcpp::Rcout << "Current int64_code: " << int64_code << std::endl;
        Rcpp::Rcout << "Last 2 bits: " << (int64_code & 3) << std::endl;
      }
      switch (int64_code & 3) {
      case 0: result.insert(result.begin(), 'A'); break;
      case 1: result.insert(result.begin(), 'C'); break;
      case 2: result.insert(result.begin(), 'T'); break;
      case 3: result.insert(result.begin(), 'G'); break;
      }
      int64_code >>= 2;
    }
    results[i] = result;
  }
  return results;
}

// [[Rcpp::export]]
Rcpp::NumericVector mutate_and_check(Rcpp::NumericVector barcode, int max_mutations, Rcpp::DataFrame r_whitelist_df, int barcode_length = 16) {
  // Convert R integer64 barcode to int64_t using RcppInt64
  std::vector<int64_t> barcodes_vec = Rcpp::fromInteger64(barcode);
  int64_t int64_barcode = barcodes_vec[0];  // Assuming you're working with a single barcode
  Rcpp::NumericVector r_whitelist = r_whitelist_df["whitelist_bcs"];
  std::unordered_set<int64_t> whitelist;
  std::vector<int64_t> whitelist_cpp = Rcpp::fromInteger64(r_whitelist);
  for (const auto& val : whitelist_cpp) {
    whitelist.insert(val);
  }
  int num_bits = 2 * barcode_length;
  std::set<int64_t> unique_results;
  recursive_flip_bits_list(int64_barcode, int64_barcode, 0, max_mutations, num_bits, unique_results, whitelist);
  // Convert the set back to a NumericVector for R using RcppInt64
  Rcpp::NumericVector results = Rcpp::toInteger64(std::vector<int64_t>(unique_results.begin(), unique_results.end()));
  return results;
}

// [[Rcpp::export]]
Rcpp::DataFrame baj_extract(std::vector<std::string>& sigstrings, 
  Rcpp::DataFrame whitelist_df, Rcpp::DataFrame df, int nthreads, 
  bool verbose = false, int max_iterations = 4, int max_mutations = 4,
  int barcode_length = 16, bool barcorrect = false, bool jaccard_on = true) {
  Rcpp::StringVector ids = df["id"];
  Rcpp::StringVector fastq_files = df["fastq_files"];
  Rcpp::StringVector qcs = df["qc"];
  std::vector<std::string> extracted_ids, extracted_barcodes, extracted_bcqcs, extracted_umis, extracted_reads, extracted_qcs;
  std::vector<std::string> updated_sig_ids;
  // Create a set for the whitelisted barcodes in bitwise representation
  Rcpp::NumericVector r_whitelist = whitelist_df["whitelist_bcs"];
  std::unordered_set<int64_t> whitelist_set;
  std::vector<int64_t> whitelist_cpp = Rcpp::fromInteger64(r_whitelist);
  for (const auto& val : whitelist_cpp) {
    whitelist_set.insert(val);
  }
  // Iterate over each signature string
  for (size_t i = 0; i < sigstrings.size(); ++i) {
    // Parse the signature string
    std::unordered_map<std::string, std::vector<int>> parsedSignature = sig_parse_ext(sigstrings[i]);
    size_t pos = sigstrings[i].find_last_of("<");
    std::string category = sigstrings[i].substr(pos);
    // Function to extract the subsequence based on parsedSignature
    auto extract_subsequence = [&](std::string id_name, std::string sequence) -> std::string {
      int start = parsedSignature[id_name][1] - 1;
      int end = parsedSignature[id_name][2];
      if (start == 0){
        start = 1;
      }
      if(end > sequence.size()){
        return category = "undecided";
      }
      if((end-start) < 5){
        return category = "undecided";
      }
      if (start >= sequence.size() || end > sequence.size()) {
        return category = "undecided";
      }
      return sequence.substr(start, end - start);
    };
    sigstrings[i] = Rcpp::as<std::string>(ids[i]) + "*" + sigstrings[i];
    if(verbose){
      Rcpp::Rcout << sigstrings[i] << std::endl;
    }
    if(category == "undecided"||category == "<undecided>"){
      continue;
    }
    if (category == "<forward>" || category == "<reverse>") {
      std::string prefix = (category == "<forward>") ? "" : "rc_";
      int start_pos = parsedSignature[prefix + "read"][1];
      int stop_pos = parsedSignature[prefix + "read"][2];
      // Set stop_pos to the length of the fastq_file string if stop_pos is zero
      if (category == "<forward>" && stop_pos == 0) {
        stop_pos = fastq_files[i].size();
      }
      // Skip this iteration if start_pos > stop_pos
      if (start_pos > stop_pos) {
        if(verbose) {
          Rcpp::Rcout << "Skipping: start_pos > stop_pos for ID " << ids[i] << std::endl;
        }
        continue;
      }
      std::string read = Rcpp::as<std::string>(fastq_files[i]).substr(start_pos, stop_pos - start_pos);
      std::string qc = Rcpp::as<std::string>(qcs[i]).substr(start_pos, stop_pos - start_pos);
      extracted_ids.push_back((sigstrings[i]) + category);
      if(category == "<reverse>"){
        std::string rc_barcode = extract_subsequence("rc_barcode", Rcpp::as<std::string>(fastq_files[i]));
        std::string revcomp_barcode = revcomp(rc_barcode);
        extracted_barcodes.push_back(revcomp_barcode);
        std::string rc_umi = extract_subsequence("rc_umi", Rcpp::as<std::string>(fastq_files[i]));
        std::string revcomp_umi = revcomp(rc_umi);
        extracted_umis.push_back(revcomp_umi);
        std::string bcqc = extract_subsequence("rc_barcode", Rcpp::as<std::string>(qcs[i]));
        extracted_bcqcs.push_back(bcqc);
      } else {
        extracted_barcodes.push_back(extract_subsequence(prefix + "barcode", Rcpp::as<std::string>(fastq_files[i])));
        extracted_umis.push_back(extract_subsequence(prefix + "umi", Rcpp::as<std::string>(fastq_files[i])));
        extracted_bcqcs.push_back(extract_subsequence(prefix + "barcode", Rcpp::as<std::string>(qcs[i])));
      }
      extracted_reads.push_back(read);
      extracted_qcs.push_back(qc);
    }
    if (category == "<rescuable_concatenates>") {
      // First split
      std::string extracted_read = extract_subsequence("read", Rcpp::as<std::string>(fastq_files[i]));
      if(extracted_read.size() < 30){
        continue;
      }
      extracted_ids.push_back(sigstrings[i] + "_splitid=<forward>");
      extracted_barcodes.push_back(extract_subsequence("barcode", Rcpp::as<std::string>(fastq_files[i])));
      extracted_umis.push_back(extract_subsequence("umi", Rcpp::as<std::string>(fastq_files[i])));
      extracted_reads.push_back(extract_subsequence("read", Rcpp::as<std::string>(fastq_files[i])));
      extracted_qcs.push_back(extract_subsequence("read", Rcpp::as<std::string>(qcs[i])));
      extracted_bcqcs.push_back(extract_subsequence("barcode", Rcpp::as<std::string>(qcs[i])));
      
      // Second split
      std::string extracted_rc_read = extract_subsequence("rc_read", Rcpp::as<std::string>(fastq_files[i]));
      if(extracted_rc_read.size() < 30){
        continue;
      }
      extracted_ids.push_back(sigstrings[i] + "_splitid=<reverse>");
      std::string rc_barcode = extract_subsequence("rc_barcode", Rcpp::as<std::string>(fastq_files[i]));
      std::string revcomp_barcode = revcomp(rc_barcode);
      extracted_barcodes.push_back(revcomp_barcode);
      std::string rc_umi = extract_subsequence("rc_umi", Rcpp::as<std::string>(fastq_files[i]));
      std::string revcomp_umi = revcomp(rc_umi);
      extracted_umis.push_back(revcomp_umi);
      extracted_reads.push_back(extract_subsequence("rc_read", Rcpp::as<std::string>(fastq_files[i])));
      extracted_qcs.push_back(extract_subsequence("rc_read", Rcpp::as<std::string>(qcs[i])));
      std::string bcqc = extract_subsequence("rc_barcode", Rcpp::as<std::string>(qcs[i]));
      extracted_bcqcs.push_back(bcqc);
    }
  }
  
  std::vector<std::string> wl_barcodes, nwl_barcodes;
  std::vector<std::string> wl_umis, nwl_umis;
  std::vector<std::string> wl_reads, nwl_reads;
  std::vector<std::string> wl_qcs, nwl_qcs;
  std::vector<std::string> wl_ids, nwl_ids;
  std::vector<std::string> wl_bcqcs, nwl_bcqcs;
  
  std::vector<int64_t> bit_barcodes = barcodes_to_bits_cpp(extracted_barcodes);
  for (size_t i = 0; i < extracted_barcodes.size(); ++i) {
    int64_t bit_barcode = bit_barcodes[i];
    if (whitelist_set.find(bit_barcode) != whitelist_set.end()) {
      wl_barcodes.push_back(extracted_barcodes[i]);
      wl_umis.push_back(extracted_umis[i]);
      wl_reads.push_back(extracted_reads[i]);
      wl_qcs.push_back(extracted_qcs[i]);
      wl_ids.push_back(extracted_ids[i]);
      wl_bcqcs.push_back(extracted_bcqcs[i]);
    } else {
      nwl_barcodes.push_back(extracted_barcodes[i]);
      nwl_umis.push_back(extracted_umis[i]);
      nwl_reads.push_back(extracted_reads[i]);
      nwl_qcs.push_back(extracted_qcs[i]);
      nwl_ids.push_back(extracted_ids[i]);
      nwl_bcqcs.push_back(extracted_bcqcs[i]);
    }
  }
  // Perform barcode error correction only if barcorrect is true
  omp_set_num_threads(nthreads);
  if (barcorrect){
    std::vector<int64_t> nwl_intcodes = barcodes_to_bits_cpp(nwl_barcodes);
    // OpenMP parallelization starts here
    #pragma omp parallel for
    for (size_t i = 0; i < nwl_barcodes.size(); ++i) {
      int64_t original_bit_barcode = nwl_intcodes[i];
      int min_hamming_distance = max_mutations;
      std::set<int64_t> min_hamming_results;
      // Ensure this function is thread-safe
      recursive_flip_bits_cpp(original_bit_barcode, original_bit_barcode, 0, max_iterations, 2 * barcode_length, whitelist_set, 
        min_hamming_distance, min_hamming_results);
      // If a single barcode is found with the minimum Hamming distance, correct the original
      if(verbose){
          #pragma omp critical
                  {
                    Rcpp::Rcout << "Original barcode: " << nwl_barcodes[i] << " Number of hamming-distance mapped barcodes found in whitelist: " << min_hamming_results.size() << std::endl;
                      if(min_hamming_results.size() <= 5){
                        std::vector<int64_t> min_hamming_results_vec(min_hamming_results.begin(), min_hamming_results.end());
                        std::vector<std::string> results = bits_to_barcodes_cpp(min_hamming_results_vec);
                        Rcpp::Rcout << "Generated barcodes: ";
                            for (const auto& str : results){
                                Rcpp::Rcout << str << ' ';
                                }
                          Rcpp::Rcout << std::endl;
    // Iterate over each generated barcode to find mismatch positions
                            for (const auto& generated_barcode : results){
                               std::vector<int> mismatch_positions = find_mismatch_positions(nwl_barcodes[i], generated_barcode);
                               Rcpp::Rcout << "Mismatch positions for " << generated_barcode << ": ";
                                  for (const auto& pos : mismatch_positions){
                                  Rcpp::Rcout << pos << ' ';
                                 }
                              Rcpp::Rcout << std::endl;
                            }
                        }
                      }
                }
      if (min_hamming_results.size() == 1) {
        int64_t corrected_bit_barcode = *min_hamming_results.begin();
        std::string corrected_barcode = bits_to_barcodes_cpp({corrected_bit_barcode})[0];
        nwl_ids[i] += "{orig_hamming_" + nwl_barcodes[i] + "_" + corrected_barcode + ":" + std::to_string(min_hamming_distance) + "}";
          if(verbose){
              #pragma omp critical
              {
                Rcpp::Rcout << "Hamming corrected barcode " << nwl_barcodes[i] << " to " << corrected_barcode << std::endl;
              }
            }
          nwl_barcodes[i] = corrected_barcode;}
      else {
        if (jaccard_on) {
        std::vector<int64_t> min_hamming_results_vec(min_hamming_results.begin(), min_hamming_results.end());
        std::vector<std::string> candidate_barcodes = bits_to_barcodes_cpp(min_hamming_results_vec);
        auto result = best_matches_jaccard_distance_qgram(nwl_barcodes[i], candidate_barcodes, 2);
        std::vector<std::string> best_matches = result.first;
        double min_jaccard_distance = result.second;
        if (best_matches.size() == 1) {
          nwl_ids[i] += "{orig_jaccard_" + nwl_barcodes[i] + "_" + best_matches[0] + ":" + std::to_string(min_jaccard_distance) + "}";
          nwl_barcodes[i] = best_matches[0];
          if (verbose) {
              #pragma omp critical
                  {
                    Rcpp::Rcout << "Jaccard corrected barcode " << nwl_barcodes[i] << " to ";
                        for (const auto& match : best_matches) {
                        Rcpp::Rcout << match << ' ';
                  }
                Rcpp::Rcout << std::endl;
                }
              }
            } else {
        if (best_matches.size() > 1) {
            nwl_ids[i] += "{jaccard_multi_" + nwl_barcodes[i] + "_" + std::to_string(best_matches.size()) + ":jac_" +
              std::to_string(min_jaccard_distance) + ":ham_" + std::to_string(min_hamming_results.size()) + "}";
              std::string concatenated_best_matches = "";
                for (const auto& match : best_matches) {
                  concatenated_best_matches += match + "|";
                }
            concatenated_best_matches.pop_back();  // Remove the trailing "|"
            nwl_barcodes[i] = concatenated_best_matches;
                }
              }
            }
          }
      }
      wl_barcodes.insert(wl_barcodes.end(), nwl_barcodes.begin(), nwl_barcodes.end());
      wl_umis.insert(wl_umis.end(), nwl_umis.begin(), nwl_umis.end());
      wl_reads.insert(wl_reads.end(), nwl_reads.begin(), nwl_reads.end());
      wl_qcs.insert(wl_qcs.end(), nwl_qcs.begin(), nwl_qcs.end());
      wl_ids.insert(wl_ids.end(), nwl_ids.begin(), nwl_ids.end());
      wl_bcqcs.insert(wl_bcqcs.end(), nwl_bcqcs.begin(), nwl_bcqcs.end());
    }
    return Rcpp::DataFrame::create(Rcpp::Named("sig_id") = wl_ids,
    Rcpp::Named("barcode") = wl_barcodes,
    Rcpp::Named("barcode_qc") = wl_bcqcs,
    Rcpp::Named("umi") = wl_umis,
    Rcpp::Named("filtered_read") = wl_reads,
    Rcpp::Named("filtered_qc") = wl_qcs);
}
