#include <Rcpp.h>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <numeric>
#include <cmath>

using namespace Rcpp;
using namespace std;
//an easy split function
std::vector<std::string> split(const std::string& s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}
// Function to parse the signature string into a data structure, returned basically [id] => [coord_vec]
std::unordered_map<std::string, std::vector<int>> sig_parse(const std::string& signature) {
  std::unordered_map<std::string, std::vector<int>> parsedSignature;
  std::stringstream ss(signature);
  std::string token;
  // first while loop breaks stuff down by pipe operator
  while (std::getline(ss, token, '|')) {
    size_t pos = token.find(":");
    std::string id = token.substr(0, pos);
    std::string coords = token.substr(pos + 1);
    std::stringstream coord_ss(coords);
    std::string coord;
    std::vector<int> coord_vec;
    // second while loop breaks stuff down by colon operator
    while (std::getline(coord_ss, coord, ':')) {
      coord_vec.push_back(std::stoi(coord));
    }
    parsedSignature[id] = coord_vec;
  }
  return parsedSignature;
}
// generate null_distance to c++-understandable
std::unordered_map<std::string, double> calculateLowerBounds(DataFrame null_distance) {
  std::unordered_map<std::string, double> lowerBounds;
  CharacterVector query_id = null_distance["query_id"];
  NumericVector null_dist = null_distance["null_distance"];
  NumericVector sd_null = null_distance["sd_null"];
  for (int i = 0; i < query_id.size(); ++i) {
    lowerBounds[as<std::string>(query_id[i])] = std::floor(null_dist[i] - std::ceil(sd_null[i]));
  } 
  return lowerBounds;
}
// generate nested unordered map keyed by [id_str][direction_str] => [exp_len][order_int][type_str]
std::unordered_map<std::string, 
  std::unordered_map<std::string,
    std::tuple<int, int, std::string>>> convertReadLayout(Rcpp::DataFrame read_layout) {
      std::unordered_map<std::string, std::unordered_map<std::string, std::tuple<int, int, std::string>>> layoutMap;
      Rcpp::StringVector id = read_layout["id"];
      Rcpp::IntegerVector expected_length = read_layout["expected_length"];
      Rcpp::StringVector type = read_layout["type"];
      Rcpp::IntegerVector order = read_layout["order"];
      Rcpp::StringVector direction = read_layout["direction"];
      for (int i = 0; i < id.size(); ++i) {
        std::string id_str = Rcpp::as<std::string>(id[i]);
        int exp_len = expected_length[i];
        std::string type_str = Rcpp::as<std::string>(type[i]);
        int order_int = order[i];
        std::string direction_str = Rcpp::as<std::string>(direction[i]);
        layoutMap[id_str][direction_str] = std::make_tuple(exp_len, order_int, type_str);
      }
      return layoutMap;
    }
//generate unordered map keyed by [id] => [order]
std::unordered_map<std::string, int> generateOrderMap(Rcpp::DataFrame read_layout) {
  std::unordered_map<std::string, int> orderMap;
  Rcpp::CharacterVector ids = read_layout["id"];
  Rcpp::IntegerVector orders = read_layout["order"];
  for (int i = 0; i < ids.size(); ++i) {
    std::string id = Rcpp::as<std::string>(ids[i]);
    int order = orders[i];
    orderMap[id] = order;
  }
  return orderMap;
}
//ID counting
std::unordered_map<std::string, int> countIDs(const std::unordered_map<std::string, std::vector<int>>& parsedSignature) {
  std::unordered_map<std::string, int> idCount;
  for (const auto& item : parsedSignature) {
    std::string id = item.first;
    const std::vector<int>& coords = item.second;
    idCount[id] += 1;
  } 
  return idCount;
} 
//adding categories to the data and filtering
std::string error_categorization(const std::string& barcoded_sigstring, 
  const std::unordered_map<std::string, double>& lowerBounds,
  const std::unordered_map<std::string, std::unordered_map<std::string, std::tuple<int, int, std::string>>>& layoutMap, 
  const std::unordered_map<std::string, int>& orderMap,
  bool verbose = false) {
  std::unordered_map<std::string, std::vector<int>> parsedSignature = sig_parse(barcoded_sigstring);
  if(verbose){
    List out;
    for (const auto& line : parsedSignature) {
      out[line.first] = line.second;
    }
    Rcpp::Rcout << "ParsedSignature output: ";
    Rf_PrintValue(out);
  }
  std::unordered_map<std::string, int> idCount = countIDs(parsedSignature);
  if(verbose){
    List out_id;
    for (const auto& line : idCount) {
      out_id[line.first] = line.second;
    }
    Rcpp::Rcout << "idCount output: ";
    Rf_PrintValue(out_id);
    Rcpp::Rcout << "Lower bounds at forw_primer: " << lowerBounds.at("forw_primer") << std::endl;
    Rcpp::Rcout << "Lower bounds at rc_forw_primer: " << lowerBounds.at("rc_forw_primer") << std::endl;
  }
  std::string category = "undecided";
  if (
      (
          // Existing conditions
          (parsedSignature.find("rc_forw_primer") == parsedSignature.end() || parsedSignature.at("rc_forw_primer")[0] >= lowerBounds.at("rc_forw_primer")) &&
            parsedSignature.at("forw_primer")[0] <= lowerBounds.at("forw_primer") &&
            parsedSignature.at("forw_primer")[0] < parsedSignature.at("rc_forw_primer")[0]
      ) ||
        ( // New condition: both primers are below the threshold
            parsedSignature.at("forw_primer")[0] <= lowerBounds.at("forw_primer") &&
              parsedSignature.at("rc_forw_primer")[0] <= lowerBounds.at("rc_forw_primer") &&
              (
                  parsedSignature.find("rev_primer") == parsedSignature.end() ||
                    parsedSignature.at("rev_primer")[0] <= lowerBounds.at("rev_primer")
              ) &&
                parsedSignature.at("forw_primer")[0] < parsedSignature.at("rev_primer")[0]
        )
  ){
    category = "forward";
    if(verbose){
      Rcpp::Rcout << "Mapped to the forward category, now attempting to find read coordinates."<< std::endl;
    }
    std::string read_key = "read";
    int read_order = orderMap.at(read_key);
    if(verbose){
      Rcpp::Rcout << "Now mapping the following sigstring: " << barcoded_sigstring << std::endl;
      Rcpp::Rcout << "Printing read_key and read order: " << read_order << " is read order & " << read_key << " is read key!" <<std::endl;
    }
    int new_start = -1;
    int new_end = -1;
    int search_distance = 1;
    std::vector<int> new_read_coords = {0, 0, 0};
    bool found_start = false, found_end = false;
    while (!found_start || !found_end) {
      if(verbose) {
        Rcpp::Rcout << "Current search distance: " << search_distance << std::endl;
      }
      for (const auto& item : parsedSignature) {
        std::string id = item.first;
        if (layoutMap.find(id) != layoutMap.end() && layoutMap.at(id).begin()->first == "forward") {
          int distance = abs(orderMap.at(id) - read_order);
          if (distance == search_distance) {
            if(verbose) {
              Rcpp::Rcout << "Checking adjacent ID: " << id << " at distance: " << distance << std::endl;
            }
            if (orderMap.at(id) < read_order) {
              if (verbose) {
                Rcpp::Rcout << "Checking id: " << id << std::endl;
              }
              if ((lowerBounds.find(id) != lowerBounds.end() && item.second[0] <= lowerBounds.at(id)) || id == "poly_a" || id == "poly_t" || id == "barcode" || id == "umi") {
                if (!found_start) {
                new_start = item.second[2] + 1;
                  }
                found_start = true; // Make sure to update the flag
                if(verbose) {
                  Rcpp::Rcout << "New start found: " << new_start << std::endl;
                }
                continue;
              }
            } else {
              if ((lowerBounds.find(id) != lowerBounds.end() && item.second[0] <= lowerBounds.at(id)) || id == "poly_a" || id == "poly_t" || id == "barcode" || id == "umi") {
                if(!found_end){
                new_end = item.second[1] - 1;
                }
                found_end = true;
                if(verbose) {
                  Rcpp::Rcout << "New end found: " << new_end << std::endl;
                }
              }
              continue;
            }
          }
        }
      }
      if (new_end == -1 && new_start != -1) {
        if(verbose){
          Rcpp::Rcout << "New end not found, setting to sentinel value 0" << std::endl;
        }
        new_end = 0;  // Note the corrected assignment operator
        break;
      }
      if (found_start && found_end) {
        if(verbose) {
          Rcpp::Rcout << "Both new start and end found, exiting loop." << std::endl;
        }
        break;
      }
      if (!found_start || !found_end) {
        search_distance++;
        if(verbose) {
          Rcpp::Rcout << "New position not found, increasing search distance to: " << search_distance << std::endl;
        }
      }
    }
    if (new_start != -1 && new_end != -1) {
      new_read_coords = {0, new_start, new_end};
      if(verbose){
        Rcpp::Rcout << "putative new coordinates are " << new_read_coords[0] << new_read_coords[1] << new_read_coords[2] << std::endl;
      }
      parsedSignature[read_key] = new_read_coords;
    }
    std::stringstream newSignature;
    for (const auto& item : parsedSignature) {
      newSignature << "|" << item.first << ":" << item.second[0] << ":" << item.second[1] << ":" << item.second[2];
    }
    newSignature << "<" << category << ">";
    if(verbose){
      Rcpp::Rcout << "putative new signature is " << newSignature.str() <<  "!" << std::endl;
    }
    return newSignature.str().substr(1);
  }
    // Check for "reverses"
    if (
        (
            // Existing conditions
            (parsedSignature.find("forw_primer") == parsedSignature.end() || parsedSignature.at("forw_primer")[0] >= lowerBounds.at("forw_primer")) &&
              parsedSignature.at("rc_forw_primer")[0] <= lowerBounds.at("rc_forw_primer") &&
              parsedSignature.at("rc_forw_primer")[0] < parsedSignature.at("forw_primer")[0]
        ) ||
          (
              // New condition: both primers are below the threshold
              parsedSignature.at("forw_primer")[0] <= lowerBounds.at("forw_primer") &&
                parsedSignature.at("rc_forw_primer")[0] <= lowerBounds.at("rc_forw_primer") &&
                (
                    parsedSignature.find("rc_rev_primer") == parsedSignature.end() ||
                      parsedSignature.at("rc_rev_primer")[0] <= lowerBounds.at("rc_rev_primer")
                ) &&
                  parsedSignature.at("rc_forw_primer")[0] < parsedSignature.at("rc_rev_primer")[0]
          )
    ){
      category = "reverse";
      if(verbose){
        Rcpp::Rcout << "Mapped to the reverse category, now attempting to find read coordinates."<< std::endl;
      }
      std::vector<int> new_read_coords = {0, 0, 0};
      std::string read_key = "rc_read";
      int read_order = orderMap.at(read_key);
      if(verbose){
        Rcpp::Rcout << "Now mapping the following sigstring: " << barcoded_sigstring << std::endl;
        Rcpp::Rcout << "Printing read_key and read order: " << read_order << " is read order & " << read_key << " is read key!" <<std::endl;
      }
      int new_start = -1;
      int new_end = -1;
      int search_distance = 1;
      bool found_start = false, found_end = false;
      while (!found_start || !found_end) {
        if(verbose) {
          Rcpp::Rcout << "Current search distance: " << search_distance << std::endl;
        }
        for (const auto& item : parsedSignature) {
          std::string id = item.first;
          if (layoutMap.find(id) != layoutMap.end() && layoutMap.at(id).begin()->first == "reverse") {
            int distance = abs(orderMap.at(id) - read_order);
            if (distance == search_distance) {
              if(verbose) {
                Rcpp::Rcout << "Checking adjacent ID: " << id << " at distance: " << distance << std::endl;
              }
              if (orderMap.at(id) < read_order) {
                if (verbose) {
                  Rcpp::Rcout << "Checking id: " << id << std::endl;
                }
                if ((lowerBounds.find(id) != lowerBounds.end() && item.second[0] <= lowerBounds.at(id)) || id == "poly_a" || id == "poly_t" || id == "barcode" || id == "umi") {
                  if (!found_end) { // Only update new_start if found_start is false
                    new_end = item.second[1] - 1;
                  }
                  found_end = true; // Make sure to update the flag
                  if(verbose) {
                    Rcpp::Rcout << "New stop found: " << new_end << std::endl;
                  }
                }
                continue;
              } else {
                if ((lowerBounds.find(id) != lowerBounds.end() && item.second[0] <= lowerBounds.at(id)) || id == "poly_a" || id == "poly_t" || id == "barcode" || id == "umi") {
                  if (!found_start) { // Only update new_start if found_start is false
                    new_start = item.second[2] + 1;
                  }
                  found_start = true;
                  if(verbose) {
                    Rcpp::Rcout << "New start found: " << new_start << std::endl;
                  }
                }
                continue;
              }
            }
          }
        }
        if (new_start == -1 && new_end != -1) {
          if(verbose){
            Rcpp::Rcout << "New start not found, setting to sentinel value 0" << std::endl;
          }
          new_start = 1;  // Note the corrected assignment operator
          break;
        }
        if (found_start && found_end) {
          if(verbose) {
            Rcpp::Rcout << "Both new start and end found, exiting loop." << std::endl;
          }
          break;
        }
        if (!found_start || !found_end) {
          search_distance++;
          if(verbose) {
            Rcpp::Rcout << "New position not found, increasing search distance to: " << search_distance << std::endl;
          }
        }
      }
      if (new_start != -1 && new_end != -1) {
        new_read_coords = {0, new_start, new_end};
        if(verbose){
          Rcpp::Rcout << "putative new coordinates are " << new_read_coords[0] <<" "<< new_read_coords[1] << " " << new_read_coords[2] << std::endl;
        }
        parsedSignature[read_key] = new_read_coords;
      }
      std::stringstream newSignature;
      for (const auto& item : parsedSignature) {
        newSignature << "|" << item.first << ":" << item.second[0] << ":" << item.second[1] << ":" << item.second[2];
      }
      newSignature << "<" << category << ">";
      if(verbose){
        Rcpp::Rcout << "putative new signature is " << newSignature.str() <<  "!" << std::endl;
      }
      return newSignature.str().substr(1);}
  // Check for "rescuable_concatenates"
   else if (parsedSignature.at("forw_primer")[0] <= lowerBounds.at("forw_primer") &&
    parsedSignature.at("rc_forw_primer")[0] <= lowerBounds.at("rc_forw_primer") && 
    parsedSignature.find("read") != parsedSignature.end() &&
    parsedSignature.find("rc_read") != parsedSignature.end()) {
    category = "rescuable_concatenates";
    if(verbose){
      Rcpp::Rcout << "Mapped as a rescuable concatenate, now attempting to rescue!" << std::endl;
      Rcpp::Rcout << barcoded_sigstring << std::endl;
    }
    // Create a vector to hold the sorted elements
    std::vector<std::tuple<std::string, int, int, int>> sortedElements;
    for (const auto& item : parsedSignature) {
      std::string id = item.first;
      if (id != "read" && id != "rc_read") {
        sortedElements.push_back(std::make_tuple(id, orderMap.at(id), item.second[1], item.second[2]));
      }
    }
    // Debug: Print elements before sorting
    if (verbose) {
      Rcpp::Rcout << "Elements before sorting:" << std::endl;
      for (const auto& item : sortedElements) {
        Rcpp::Rcout << std::get<0>(item) << " : " << std::get<2>(item) << "-" << std::get<3>(item) << std::endl;
      }
    }
    // Sort by start position
    std::sort(sortedElements.begin(), sortedElements.end(),
      [](const auto& a, const auto& b) {
        return std::get<2>(a) < std::get<2>(b);
      });
    // Debug: Print elements after sorting
    if (verbose) {
      Rcpp::Rcout << "Elements after sorting:" << std::endl;
      for (const auto& item : sortedElements) {
        Rcpp::Rcout << std::get<0>(item) << " : " << std::get<2>(item) << "-" << std::get<3>(item) << std::endl;
      }
    }
    int read_start = parsedSignature.at("read")[2] + 1;
    int read_end = parsedSignature.at("rc_read")[1] - 1;
    if(verbose){
      Rcpp::Rcout << "read start is " << read_start << " and read end is " << read_end << std::endl;
    }
    // Find read start and end based on sorted elements
    for (const auto& item : sortedElements) {
      std::string id = std::get<0>(item);
      int start = std::get<2>(item);
      int end = std::get<3>(item);
      if (start > read_start) {
        read_end = start - 1;
        break;
      }
      if (end < read_end) {
        read_start = end + 1;
      }
    }
    // Debug: Print new read start and end
    if (verbose) {
      Rcpp::Rcout << "New read start: " << read_start << ", New read end: " << read_end << std::endl;
    }
    // Update parsedSignature for read and rc_read
    parsedSignature["read"] = {0, read_start, read_end};
    parsedSignature["rc_read"] = {0, read_end + 1, parsedSignature.at("rc_read")[2]};
    // Step 4: Re-sort the complete list, including updated read and rc_read
    sortedElements.push_back(std::make_tuple("read", orderMap.at("read"), read_start, read_end));
    sortedElements.push_back(std::make_tuple("rc_read", orderMap.at("rc_read"), read_end + 1, parsedSignature.at("rc_read")[2]));
    std::sort(sortedElements.begin(), sortedElements.end(),
      [](const auto& a, const auto& b) {
        return std::get<2>(a) < std::get<2>(b);
      });
    if (verbose) {
      Rcpp::Rcout << "Elements after sorting:" << std::endl;
      for (const auto& item : sortedElements) {
        Rcpp::Rcout << std::get<0>(item) << " : " << std::get<2>(item) << "-" << std::get<3>(item) << std::endl;
      }
    }
    // Initialize variables for the new start and stop positions
    int new_read_start = -1;
    int new_read_end = -1;
    int new_rc_read_start = -1;
    int new_rc_read_end = -1;
    // Iterate through the sorted elements to update start and stop positions
    for (size_t i = 0; i < sortedElements.size(); ++i) {
      if (std::get<0>(sortedElements[i]) == "read") {
        new_read_start = std::get<3>(sortedElements[i - 1]) + 1;  // Stop of the element before 'read' + 1
        new_read_end = std::get<2>(sortedElements[i + 1]) - 1;  // Start of the element after 'read' - 1
      }
      if (std::get<0>(sortedElements[i]) == "rc_read") {
        new_rc_read_start = std::get<3>(sortedElements[i - 1]) + 1;  // Stop of the element before 'rc_read' + 1
        new_rc_read_end = std::get<2>(sortedElements[i + 1]) - 1;  // Start of the element after 'rc_read' - 1
      }
    }
    // Update parsedSignature with new start and stop positions
    parsedSignature["read"] = {0, new_read_start, new_read_end};
    parsedSignature["rc_read"] = {0, new_rc_read_start, new_rc_read_end};
    // Debug: Print new read and rc_read start and end
    if (verbose) {
      Rcpp::Rcout << "New read start: " << new_read_start << ", New read end: " << new_read_end << std::endl;
      Rcpp::Rcout << "New rc_read start: " << new_rc_read_start << ", New rc_read end: " << new_rc_read_end << std::endl;
    }
    // Find and update 'read' and 'rc_read' in sortedElements
    for (auto& item : sortedElements) {
      if (std::get<0>(item) == "read") {
        std::get<2>(item) = new_read_start;
        std::get<3>(item) = new_read_end;
      }
      if (std::get<0>(item) == "rc_read") {
        std::get<2>(item) = new_rc_read_start;
        std::get<3>(item) = new_rc_read_end;
      }
    }
    // Generate the new concatenated signature
    std::stringstream newSignature;
    for (const auto& item : sortedElements) {
      newSignature << "|" << std::get<0>(item) << ":" << 0 << ":" << std::get<2>(item) << ":" << std::get<3>(item);
    }
    // Debug: Print new concatenated signature
    if (verbose) {
      Rcpp::Rcout << "New concatenated signature: " << newSignature.str() << std::endl;
    }
    newSignature << "<" << category << ">";
    return newSignature.str().substr(1);
  }
  // Check for "parallel_f_concatenates"
  if (category == "forwards" && idCount["forw_primer"] > 1) {
    category = "parallel_f_concatenates";
  }
  // Check for "parallel_r_concatenates"
  if (category == "reverses" && idCount["rc_forw_primer"] > 1) {
    category = "parallel_r_concatenates";
  }
  return barcoded_sigstring + "<" + category + ">";
}
//findRead statically annotates read in the string itself and sorts everything in the read sequentially
std::string read_annotator(const std::string& signatureString, const std::unordered_map<std::string, int>& orderMap) {
  std::vector<std::pair<int, std::string>> orderedComponents;
  // Split the signature string by '|'
  std::vector<std::string> components = split(signatureString, '|');
  for (const auto& component : components) {
    // Extract the ID from each component (the part before the first ':')
    std::string id = component.substr(0, component.find(":"));
    // Check if this ID is in our order map
    if (orderMap.find(id) != orderMap.end()) {
      int order = orderMap.at(id);
      orderedComponents.push_back({order, component});
    }
  }
  // Sort by the order
  std::sort(orderedComponents.begin(), orderedComponents.end());
  std::stringstream reorderedSigString;
  int order_of_read = orderMap.at("read");
  int order_of_rc_read = orderMap.at("rc_read");
  for (const auto& item : orderedComponents) {
    reorderedSigString << "|" << item.second;
    // If the order of the current item is one less than the order of "read", then insert "read"
    if (item.first == order_of_read-1) {
      // Inserting a read with arbitrary start and stop positions (0, 100) for now
      // 
      reorderedSigString << "|read:0:0:100";
    }
    // If the order of the current item is one less than the order of "rc_read", then insert "rc_read"
    if (item.first == order_of_rc_read-1) {
      // Inserting an rc_read with arbitrary start and stop positions (0, 100) for now
      reorderedSigString << "|rc_read:0:0:100";
    }
  }
  return reorderedSigString.str().substr(1);  // Remove the leading "|"
}
//read batching & error mapping
std::string bajbatch(const std::string sigstring, 
  const std::unordered_map<std::string, std::unordered_map<std::string, std::tuple<int, int, std::string>>>& layoutMap, 
  const std::unordered_map<std::string, int>& orderMap,
  const std::unordered_map<std::string, double>& lowerBounds, 
  bool verbose = false) {
  std::stringstream barcoded_sigstring;
  std::unordered_map<std::string, std::vector<int>> parsed_sigstring = sig_parse(sigstring);
  for (const auto& annot : parsed_sigstring) {
    std::string id = annot.first;
    std::vector<int> coords = annot.second;
    int edit_distance = coords[0];
    int start_pos = coords[1];
    int stop_pos = coords[2];
    if(verbose) {
      Rcpp::Rcout << "Checking piece identity!"<< std::endl;
      Rcpp::Rcout << id << std::endl;
    }
    // Skip this iteration if the edit distance is greater than the lower bound
    if (lowerBounds.find(id) != lowerBounds.end() && edit_distance > lowerBounds.at(id)) {
      barcoded_sigstring << "|" << id << ":" << edit_distance << ":" << start_pos << ":" << stop_pos;
      if(verbose) {
        Rcpp::Rcout << "Skipping this id: " << id << std::endl;
        Rcpp::Rcout << "Updated Signature so far: " << barcoded_sigstring.str() << std::endl;
      }
      continue;
    }
    // Check if this ID is a flanking adapter
    if (layoutMap.find(id) != layoutMap.end()) {
      std::string direction = layoutMap.at(id).begin()->first;
      std::string type = std::get<2>(layoutMap.at(id).begin()->second);
      if (type == "flanking_adapter") {
        std::string barcode_id = (direction == "forward") ? "barcode" : "rc_barcode";
        std::string umi_id = (direction == "forward") ? "umi" : "rc_umi";
        // Annotate barcode
        if (layoutMap.find(barcode_id) != layoutMap.end() && layoutMap.at(barcode_id).find(direction) != layoutMap.at(barcode_id).end()) {
          int barcode_length = std::get<0>(layoutMap.at(barcode_id).at(direction));
          int new_start, new_stop;
          if (direction == "forward") {
            new_start = stop_pos + 1;
            new_stop = new_start + barcode_length - 1;
          }
          if (direction == "reverse") {
            new_stop = start_pos - 1;
            new_start = new_stop - barcode_length + 1;
          }
          barcoded_sigstring << "|" << barcode_id << ":" << 0 << ":" << new_start << ":" << new_stop;
        }
        // Annotate UMI
        if (layoutMap.find(umi_id) != layoutMap.end() && layoutMap.at(umi_id).find(direction) != layoutMap.at(umi_id).end()) {
          int barcode_length = std::get<0>(layoutMap.at(barcode_id).at(direction));
          int umi_length = std::get<0>(layoutMap.at(umi_id).at(direction));
          int new_start, new_stop;
          
          if (direction == "forward") {
            new_start = stop_pos + barcode_length + 1;
            new_stop = new_start + umi_length - 1;
          } else {
            new_stop = start_pos - barcode_length - 1;
            new_start = new_stop - umi_length + 1;
          }
          barcoded_sigstring << "|" << umi_id << ":" << 0 << ":" << new_start << ":" << new_stop;
        }
      } 
    }
    // Append the existing id and coordinates
    barcoded_sigstring << "|" << id << ":" << edit_distance << ":" << start_pos << ":" << stop_pos;
    if(verbose){
      Rcpp::Rcout << "Updated Signature with barcodes and umis:  " << barcoded_sigstring.str().substr(1) << std::endl;
    }
  }
  std::string ordered_sigstring = read_annotator(barcoded_sigstring.str().substr(1), orderMap);
  if(verbose){
    std::string category = error_categorization(ordered_sigstring, lowerBounds,layoutMap, orderMap, true);
    return category;
  } else {
    std::string category = error_categorization(ordered_sigstring, lowerBounds,layoutMap, orderMap, false);
    return category;
  }
}

// [[Rcpp::export]]
StringVector bajbatch(DataFrame null_distance, DataFrame read_layout, StringVector sigstrings, bool verbose = false) {
  if (verbose) {
    Rcpp::Rcout << "Null Distance DataFrame: " << std::endl;
    CharacterVector query_id = null_distance["query_id"];
    NumericVector null_dist = null_distance["null_distance"];
    NumericVector sd_null = null_distance["sd_null"];
    for (int i = 0; i < query_id.size(); ++i) {
      Rcpp::Rcout << as<std::string>(query_id[i]) << ", " << null_dist[i] << ", " << sd_null[i] << std::endl;
    }
  }
  // Calculate the lower bounds for each id
  std::unordered_map<std::string, double> lowerBounds = calculateLowerBounds(null_distance);
  if(verbose) {
    Rcpp::Rcout << "lowerBounds: " << std::endl;
    for (const auto& pair : lowerBounds) {
      Rcpp::Rcout << "ID: " << pair.first 
                  << ", Value: " << pair.second << std::endl;
    }
  }
  // Convert read_layout to an appropriate C++ data structure
  std::unordered_map<std::string, std::unordered_map<std::string, 
    std::tuple<int, int, std::string>>> layoutMap = convertReadLayout(read_layout);  
  std::unordered_map<std::string, int> orderMap = generateOrderMap(read_layout);
  if (verbose) {
    Rcpp::Rcout << "Layout Map: " << std::endl;
    for (const auto& outer_pair : layoutMap) {
      for (const auto& inner_pair : outer_pair.second) {
        Rcpp::Rcout << "ID: " << outer_pair.first 
                    << ", Direction: " << inner_pair.first
                    << ", Expected Length: " << std::get<0>(inner_pair.second)
                    << ", Order: " << std::get<1>(inner_pair.second) << std::endl
                    << ", Type: " << std::get<2>(inner_pair.second) << std::endl;
      }
    }
  }
  // Iterate through each signature string to parse and update it
  for (int i = 0; i < sigstrings.size(); ++i) {
    std::string sigstring = as<std::string>(sigstrings[i]);
    if(verbose){
      std::string output = bajbatch(sigstring, layoutMap, orderMap, lowerBounds, true);
      sigstrings[i] = output;
    } else {
      std::string output = bajbatch(sigstring, layoutMap, orderMap, lowerBounds, false);
      sigstrings[i] = output;
    }
  }
  return sigstrings;
}
