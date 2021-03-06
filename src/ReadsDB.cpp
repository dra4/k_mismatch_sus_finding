#include "ReadsDB.hpp"
#include "util.hpp"

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <climits>

#include <dirent.h>
#include <assert.h>

// UTILITY FUNCTIONS -------------------------------------------------

int getFilenames(const std::string& directory,
                 std::vector< std::string >& fileNames) {
    DIR *dir;
    struct dirent *ent;
    dir = opendir(directory.c_str());
    if (dir != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            std::string ename = ent->d_name;
            if(ename.size() > 3 &&
               ename.substr(ename.size()-3) == ".fa") {
                // printf("%s\n", ent->d_name);
                fileNames.push_back(directory + "/" + ent->d_name);
            }
            if(ename.size() > 3 &&
               ename.substr(ename.size()-3) == ".fq") {
                // printf("%s\n", ent->d_name);
                fileNames.push_back(directory + "/" + ent->d_name);
            }
        }
        closedir(dir);
    }
    else {
        std::cout << "Failed to read directory " << EXIT_FAILURE
                  << std::endl;
        return 1;
    }
    return 0;
}


// READS DATABASE ---------------------------------------------------

// clip the read and add it to the database
void ReadsDB::addRead(std::string& read, std::string& header, preads_t& freads,
                      rlengths_t& flengths, rnames_t& fnames) {
    if(read.length() == 0)
        return;

    // skipping first x character and last y characters
    freads.push_back(read);
    fnames.push_back(header);
    flengths.push_back(read.length());

}

// read fasta file and add it to the read set; update lengths too
void ReadsDB::readFasta(std::string& fname, preads_t& freads,
                        rlengths_t& flengths, rnames_t& fnames) {
    std::ifstream fs(fname);
    std::string line = "";
    std::string read = "";
    std::string header = "";

    while(getline(fs, line)) {
        if(line[0] == '>') {
            addRead(read, header, freads, flengths, fnames);
            read = "";
            line = line.substr(1);
            header = trim(line);
        } else if(is_valid_apha(line[0])) {
            line  = trim(line);
            read += line;
        }
    }
    addRead(read, header, freads, flengths, fnames);
    assert(freads.size() == flengths.size());
    assert(freads.size() == fnames.size());
}

// read fasta file and add it to the read set; update lengths too
void ReadsDB::readFastq(std::string& fname, preads_t& freads,
                        rlengths_t& flengths, rnames_t& fnames) {
    std::ifstream fs(fname);
    std::string line = "";
    std::string read = "";
    std::string header = "";

    while(getline(fs, line)) {
        if(line[0] == '@') {
            read = "";
            line = line.substr(1);
            header = trim(line);
        } else if(line[0] == 'A' || line[0] == 'C' ||
                  line[0] == 'G' || line[0] == 'T') {
            line  = trim(line);
            read += line;
        } else if(line[0] == '+') {
            addRead(read, header, freads, flengths, fnames);
            read = "";
            header = "";
        }
    }
    // addRead(read, header, freads, flengths, fnames);
    assert(freads.size() == flengths.size());
    assert(freads.size() == fnames.size());
}


// find the maximum length of sequence in the read-set
unsigned ReadsDB::findMax(){
    unsigned readMaxLen = 0;
    for(unsigned int ix = 0; ix < m_readsStore.size(); ix++){
        for(unsigned int iy = 0; iy < m_readsStore[ix].size(); iy++){
            if(readMaxLen < m_readsStore[ix][iy].length())
                readMaxLen = m_readsStore[ix][iy].length();
        }
    }
    return readMaxLen;
}

// find the minimum length of sequence in the read-set
unsigned ReadsDB::findMin(){
    unsigned readMinLen = UINT_MAX;
    for(unsigned int ix = 0; ix < m_readsStore.size(); ix++){
        for(unsigned int iy = 0; iy < m_readsStore[ix].size(); iy++){
            if(readMinLen > m_readsStore[ix][iy].length())
                readMinLen = m_readsStore[ix][iy].length();
        }
    }
    return readMinLen;
}

// pad pChar to every read in the read set,
//   such that all reads have length mxLength + 1
void ReadsDB::padStrings(const unsigned& mxLength, const char& pChar){
    for(unsigned int ix = 0; ix < m_readsStore.size(); ix++){
        for(unsigned int iy = 0; iy < m_readsStore[ix].size(); iy++){
            unsigned fillIn = mxLength - m_readsStore[ix][iy].length();
            for(unsigned k = 0; k < fillIn; k++)
                m_readsStore[ix][iy] += pChar;
        }
    }
}

// pad reads with fillChar to all reads in the db,
//    such that all reads have length m_readMaxLength + 1
unsigned ReadsDB::padReads(const char& fillChar){
    m_readMaxLength = findMax();

#ifdef DEBUG_RDB
    // print ONLY for debugging!
    std::cout << "\"NMAX_LEN\" : " << m_readMaxLength << "," << std::endl
              << "\"NMIN_LEN\" : " << findMin() << "," << std::endl;
#endif
    padStrings(m_readMaxLength + 1, fillChar);

#ifdef DEBUG_RDB
    // print ONLY for debugging!
    std::cout << "\"PMAX_LEN\" : " << findMax() << "," << std::endl
              << "\"PMIN_LEN\" : " << findMin() << "," << std::endl;
#endif
    // m_readMaxLength += 1;
    return m_readMaxLength;
}

// return the read database as a single string
void ReadsDB::asString(std::string& total){
    unsigned tx = 0;
    assert(m_nReads > 0);
    assert(m_readMaxLength > 0);
    total.resize(m_nReads * (m_readMaxLength + 1));

    for(unsigned ix = 0; ix < m_readsStore.size(); ix++){
        for(unsigned iy = 0; iy < m_readsStore[ix].size(); iy++){
            for(unsigned iz = 0; iz < m_readsStore[ix][iy].length(); iz++){
                total[tx] = m_readsStore[ix][iy][iz];
                tx++;
            }
        }
    }
}

// convert read database as a single string m_dbString
void ReadsDB::buildDBS(){
    // TODO: clear the read store after building the db string
    asString(m_dbString);
}

// For debugging and output
void ReadsDB::writeMeta(std::ostream& ots){
    ots << "\"nfiles\"      : " << m_readFiles.size() << "," << std::endl;

    ots << "\"total_reads\" : " << m_nReads
        << "," << std::endl;
#ifdef DEBUG_RDB
    ots << "\"n_read_counts\" : " << m_readCounts.size()
        << "," << std::endl;
#endif
    ots << "\"file_names\" : [" << std::endl;
    for(unsigned i = 0; i < m_readFiles.size(); i++){
        ots << " \t\"" << m_readFiles[i] << "\"," << std::endl;
    }
    ots << " \t\"\"" << std::endl << " ]," << std::endl;

    ots << "\"read_counts\" : [" << std::endl;
    for (unsigned i = 0; i < m_readCounts.size(); i++) {
        ots << " \t[" << i << "," << m_readCounts[i]
#ifdef DEBUG_RDB
            << "," << m_readLengths[i].size()
            << "," <<  m_readCtsPfxSum[i]
#endif
            << "]," << std::endl;
    }
    ots << " \t[]" << std::endl << " ]," << std::endl;
}

void ReadsDB::writeNames(std::ostream& ots){
    ots << "\"read_names\" : [" << std::endl;
    unsigned idx = 0;
    for(unsigned i = 0; i < m_readNames.size(); i++){
        for(unsigned j = 0; j < m_readNames[i].size(); j+=2) {
            // ots << " \t\"" << idx << "\": \"" << m_readNames[i][j]
            ots << " \t\"" << m_readNames[i][j]
                << "\"," << std::endl;
            idx++;
        }
    }
    // ots << " \t\"" << idx << "\": \"\"" << " }," << std::endl;
    ots << "\t \"\" ]," << std::endl;
}

void ReadsDB::init(char padChar){
    m_nReads = 0;
    std::sort(m_readFiles.begin(), m_readFiles.end());
    m_readsStore.resize(m_readFiles.size());
    m_readCounts.resize(m_readFiles.size());
    m_readCtsPfxSum.resize(m_readFiles.size());
    m_readLengths.resize(m_readFiles.size());
    m_readNames.resize(m_readFiles.size());

    for(unsigned i  = 0; i <  m_readFiles.size(); i++){
        if(ends_with(m_readFiles[i], ".fa") ||
           ends_with(m_readFiles[i], ".tfa") ||
           ends_with(m_readFiles[i], ".fas") ||
           ends_with(m_readFiles[i], ".fasta")) {
            readFasta(m_readFiles[i], m_readsStore[i],
                      m_readLengths[i], m_readNames[i]);
        }
        else if(ends_with(m_readFiles[i], ".fq")){
            readFastq(m_readFiles[i], m_readsStore[i],
                      m_readLengths[i], m_readNames[i]);
        }
        m_readCtsPfxSum[i] = m_readCounts[i] = m_readsStore[i].size();
        m_nReads += m_readsStore[i].size();
    }

    // Pad the strings upto maximum length
    if(padChar != 0)
        m_readMaxLength = padReads(padChar);
    else
        m_readMaxLength = findMax();

    //  Prefix sum read counts
    for(unsigned i = 1; i < m_readCtsPfxSum.size(); i++)
        m_readCtsPfxSum[i] += m_readCtsPfxSum[i-1];

}

// initilaize reads database with a list of input files
ReadsDB::ReadsDB(std::vector<std::string>& fnames, char padChar):m_readFiles(fnames){
    init(padChar);
    // writeMeta(std::cout);
}

ReadsDB::ReadsDB(std::string& dirname, char padChar){
    // Construct the Reads Store
    getFilenames(dirname, m_readFiles);
    init(padChar);
    // writeMeta(std::cout);
}

// Initialize the reads database using the config params
ReadsDB::ReadsDB(std::string& dirname, std::vector<std::string>& fnames,
                 char padChar){
    if(dirname.length() > 0){
        getFilenames(dirname, m_readFiles);
    }
    if(fnames.size() > 0) {
        for(auto it = fnames.begin();it != fnames.end();++it)
            m_readFiles.push_back(*it);
    }
    init(padChar);
    //writeMeta(std::cout);
}

// Lookup the file ID corresponding to read ID
int ReadsDB::getFileOfRead(const unsigned& readID) const{
    // TODO: change to binary search
    for(int i = 0; i < (int)m_readCtsPfxSum.size(); i++)
        if(m_readCtsPfxSum[i] > readID) // TODO:: should this be >=?
            return i;
    return -1;
}

// Lookup the file ID correspond to DB String position
unsigned ReadsDB::getDBSFileId(const unsigned& dbStrPosition){
    unsigned readID = dbStrPosition / (m_readMaxLength + 1);
    int fid =  getFileOfRead(readID);
    assert(fid >= 0);
    return (unsigned)fid;
}

// Lookup the read ID correspond to DB String position
unsigned ReadsDB::getDBSReadId(const unsigned& dbStrPosition){
    // TODO: handle this appropriately
    unsigned readID = dbStrPosition / (m_readMaxLength + 1);
    int fofRead = getFileOfRead(readID);
    assert(fofRead >= 0);
    assert(fofRead < ((int)m_readFiles.size()));
    if(fofRead > 0){
        if((readID - m_readCtsPfxSum[fofRead - 1]) >=
           m_readLengths[fofRead].size()) {
            std::cout << fofRead  << " " << readID << " "
                      << dbStrPosition << " " << m_readMaxLength << " "
                      << m_readLengths[fofRead].size() << std::endl;
            std::flush(std::cout);
        }
        assert(readID >= m_readCtsPfxSum[fofRead - 1]);
        readID -=  m_readCtsPfxSum[fofRead - 1];
        assert(readID < m_readLengths[fofRead].size());
    }
    return readID;
}

// Lookup the read Pos correspond to DB String position
unsigned ReadsDB::getDBSReadPos(const unsigned& dbStrPosition){
    return (dbStrPosition) % (m_readMaxLength + 1);
}

// Lookup the char before the given DB String position
char ReadsDB::getDBSCharBefore(const unsigned& dbStrPosition){
    assert(dbStrPosition < m_dbString.size());

    if(dbStrPosition == 0)
        return '$';
    return m_dbString[dbStrPosition - 1];
}

// Lookup the char at the given DB String position
char ReadsDB::getDBSCharAt(const unsigned& dbStrPosition){
    assert(dbStrPosition < m_dbString.size());

    return m_dbString[dbStrPosition];
}

// Lookup the read length, given fileID and readID
unsigned ReadsDB::getReadLength(const unsigned& fileID, const unsigned& readID){
    assert(fileID < m_readFiles.size());
    if(readID > m_readLengths[fileID].size()){
        std::cout << "FID " << fileID << " RID " << readID
                  << " RLEN "<< m_readLengths[fileID].size()
                  << std::endl;
        std::flush(std::cout);
    }
    assert(readID < m_readLengths[fileID].size());
    return m_readLengths[fileID][readID];
}

const std::string& ReadsDB::getRead(const unsigned& fileID, const unsigned& readID){
    assert(fileID < m_readFiles.size());
    if(readID > m_readLengths[fileID].size()){
        std::cout << "FID " << fileID << " RID " << readID
                  << " RLEN "<< m_readLengths[fileID].size()
                  << std::endl;
        std::flush(std::cout);
    }
    assert(readID < m_readLengths[fileID].size());
    return m_readsStore[fileID][readID];
}

unsigned ReadsDB::getReadId(const unsigned& fileID, const unsigned& readID){
    assert(fileID < m_readFiles.size());
    assert(readID < m_readLengths[fileID].size());
    return (fileID == 0 ? 0 : m_readCtsPfxSum[fileID - 1]) + readID;
}

const std::string& ReadsDB::getReadById(const unsigned& grID) const{
    int fofRead = getFileOfRead(grID);
    assert(fofRead >= 0);
    assert(fofRead < ((int)m_readFiles.size()));
    unsigned readID = grID;
    if(fofRead > 0)
        readID -= m_readCtsPfxSum[fofRead - 1];
    // std::cout << grID << std::endl;
    assert(readID >= 0);
    assert(readID < m_readsStore[fofRead].size());
    return m_readsStore[fofRead][readID];
}

const std::string& ReadsDB::getReadNameById(const unsigned& grID) const{
    int fofRead = getFileOfRead(grID);
    assert(fofRead >= 0);
    assert(fofRead < ((int)m_readFiles.size()));
    unsigned readID = grID;
    if(fofRead > 0)
        readID -= m_readCtsPfxSum[fofRead - 1];
    assert(readID >= 0);
    assert(readID < m_readNames[fofRead].size());
    return m_readNames[fofRead][readID];
}

unsigned ReadsDB::getReadLengthById(const unsigned& grID){
    int fofRead = getFileOfRead(grID);
    assert(fofRead >= 0);
    assert(fofRead < ((int)m_readFiles.size()));
    unsigned readID = grID;
    if(fofRead > 0)
        readID -= m_readCtsPfxSum[fofRead - 1];
    assert(readID >= 0);
    assert(readID < m_readLengths[fofRead].size());
    return m_readLengths[fofRead][readID];
}

const std::string& ReadsDB::getFileName(const unsigned& idx){
    assert(idx < m_readFiles.size());
    return m_readFiles[idx];
}
