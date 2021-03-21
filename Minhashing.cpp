#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <string.h>
#include <math.h>
#include <vector>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include<bits/stdc++.h>
#include <utility>
#include <type_traits>
#include <typeinfo>
#include <cxxabi.h>
#include <chrono>
#include "Shingles.hpp"
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <unordered_map>
#include "characteristic_matrix.h"
#include "minhash_signature.h"
using namespace std;
using namespace std::chrono;

const int k = 20;  ///{k-shingle size: 20, 30, 40}.
int doc_cont = 0;
const int signature_size = 36;
const int bandas = 6;
const int filas_por_banda = signature_size/bandas;
const double t = 0.75; //threshold
typedef unsigned int Signature;

vector<string> k_shingles;
Signature minhashing;
vector<Signature> vector_signature;
vector<vector<Signature>> signature_matrix;
vector<Signature> signature_band;
map<int, set<int>> buckets_banda_0;
map<int, set<int>> buckets_banda_1;
map<int, set<int>> buckets_banda_2;
map<int, set<int>> buckets_banda_3;
map<int, set<int>> buckets_banda_4;
map<int, set<int>> buckets_banda_5;
map<pair<int, int>, double> candidatos;

////////////////////////////////////////////////////////

Shingles::Shingles(){
    A = (rand() % bigPrimeNumber) + 1;
    B = (rand() % bigPrimeNumber) + 1;
}
Shingles::~Shingles(){}

Shingles::Signature

Shingles::compute(const std::vector<std::string> &seqs)const{
    std::hash<std::string> stringHash;
    Signature minShingleHash = bigPrimeNumber;

    for (auto &val : seqs){
        std::size_t hashedWord = stringHash(val);
        std::size_t shingleID = hashedWord;
        Signature shingleHash = (((unsigned long)A * (unsigned long)shingleID) + B) % bigPrimeNumber;
        if (shingleHash < minShingleHash){
            minShingleHash = shingleHash;
        }
    }

    return minShingleHash;
}

////////////////////////////////////////////////////////

double jaccardSimilarity(int d1, int d2){
    double value = 0;
    for (int k = 0; k < signature_size; k++){
        if (signature_matrix[d1][k] == signature_matrix[d2][k]){
            value++;
        }
    }
    value = value / signature_size;
    return value;
}

/////////////////////////////////////////////////////////

void generarCandidatos(){
    cout << "|||Se generan los candidatos mediante:" << endl;
    hash<string> h;
    Signature aux;
    string aux3;
    int cont_pares_candidatos = 0;
    cout << "Hashing de bandas atraves de las columnas..." << endl;
    for (int banda = 0; banda < bandas; banda++){
        for (int doc = 0; doc < (doc_cont-1); doc++){
            for (int row = 0; row < filas_por_banda; row++){
                aux = signature_matrix[doc][banda * filas_por_banda + row];
                string aux2 = to_string(aux);
                aux3 = aux3 + aux2;
            }
            int hashed = h(aux3);
            aux3.clear();
            if(banda == 0){
                buckets_banda_0[hashed].insert(doc);
            }
            else if(banda == 1){
                buckets_banda_1[hashed].insert(doc);
            }
            else if(banda == 2){
                buckets_banda_2[hashed].insert(doc);
            }
            else if(banda == 3){
                buckets_banda_3[hashed].insert(doc);
            }
            else if(banda == 4){
                buckets_banda_4[hashed].insert(doc);
            }
            else{
                buckets_banda_5[hashed].insert(doc);
            }
        }
        if(banda==0){
           /*  cout << buckets_banda_0.size() << endl; */
            map<int, set<int>>::iterator it;
            set<int>::iterator it2;
            set<int>::iterator it3;
            for (it = buckets_banda_0.begin(); it != buckets_banda_0.end(); ++it){
                if ((it->second).size() > 1){
                    for (auto it2 = (it->second).begin(), e = (it->second).end(); it2 != e; it2++){
                        cout << "" << *it2 << " ";
                        for (auto it3 = next(it2); it3 != e; it3++){
                            pair<int, int> p;
                            p.first = *it2;
                            p.second = *it3;
                            candidatos[p] = 0;
                            cont_pares_candidatos++;
                        }
                    }
                }
            }

            cout << "\n\n" << "Numero de candidatos en banda 0: " << cont_pares_candidatos << endl;
            cont_pares_candidatos = 0;
        }
        else if(banda==1){
            /* cout << buckets_banda_1.size() << endl; */
            map<int, set<int>>::iterator it;
            set<int>::iterator it2;
            set<int>::iterator it3;
            for (it = buckets_banda_1.begin(); it != buckets_banda_1.end(); ++it){
                if ((it->second).size() > 1){
                    for (auto it2 = (it->second).begin(), e = (it->second).end(); it2 != e; it2++){
                        cout << "" << *it2 << " ";
                        for (auto it3 = next(it2); it3 != e; it3++){
                            pair<int, int> p;
                            p.first = *it2;
                            p.second = *it3;
                            candidatos[p] = 0;
                            cont_pares_candidatos++;
                        }
                    }
                }
            }
            cout << "\n\n" <<"Numero de candidatos en banda 1: " << cont_pares_candidatos << endl;
            cont_pares_candidatos = 0;
        }
        else if(banda==2){
            /* cout << buckets_banda_1.size() << endl; */
            map<int, set<int>>::iterator it;
            set<int>::iterator it2;
            set<int>::iterator it3;
            for (it = buckets_banda_2.begin(); it != buckets_banda_2.end(); ++it){
                if ((it->second).size() > 1){
                    for (auto it2 = (it->second).begin(), e = (it->second).end(); it2 != e; it2++){
                        cout << "" << *it2 << " ";
                        for (auto it3 = next(it2); it3 != e; it3++){
                            pair<int, int> p;
                            p.first = *it2;
                            p.second = *it3;
                            candidatos[p] = 0;
                            cont_pares_candidatos++;
                        }
                    }
                }
            }
            cout << "\n\n" <<"Numero de candidatos en banda 1: " << cont_pares_candidatos << endl;
            cont_pares_candidatos = 0;
        }
        else if(banda==3){
            /* cout << buckets_banda_1.size() << endl; */
            map<int, set<int>>::iterator it;
            set<int>::iterator it2;
            set<int>::iterator it3;
            for (it = buckets_banda_3.begin(); it != buckets_banda_3.end(); ++it){
                if ((it->second).size() > 1){
                    for (auto it2 = (it->second).begin(), e = (it->second).end(); it2 != e; it2++){
                        cout << "" << *it2 << " ";
                        for (auto it3 = next(it2); it3 != e; it3++){
                            pair<int, int> p;
                            p.first = *it2;
                            p.second = *it3;
                            candidatos[p] = 0;
                            cont_pares_candidatos++;
                        }
                    }
                }
            }
            cout << "\n\n" <<"Numero de candidatos en banda 1: " << cont_pares_candidatos << endl;
            cont_pares_candidatos = 0;
        }
        else if(banda==4){
            /* cout << buckets_banda_1.size() << endl; */
            map<int, set<int>>::iterator it;
            set<int>::iterator it2;
            set<int>::iterator it3;
            for (it = buckets_banda_4.begin(); it != buckets_banda_4.end(); ++it){
                if ((it->second).size() > 1){
                    for (auto it2 = (it->second).begin(), e = (it->second).end(); it2 != e; it2++){
                        cout << "" << *it2 << " ";
                        for (auto it3 = next(it2); it3 != e; it3++){
                            pair<int, int> p;
                            p.first = *it2;
                            p.second = *it3;
                            candidatos[p] = 0;
                            cont_pares_candidatos++;
                        }
                    }
                }
            }
            cout << "\n\n" <<"Numero de candidatos en banda 1: " << cont_pares_candidatos << endl;
            cont_pares_candidatos = 0;
        }
        else{
            /* cout << buckets_banda_2.size() << endl; */
            map<int, set<int>>::iterator it;
            set<int>::iterator it2;
            set<int>::iterator it3;
            for(it = buckets_banda_5.begin(); it != buckets_banda_5.end(); ++it){
                if((it->second).size() > 1){
                    for (auto it2 = (it->second).begin(), e = (it->second).end(); it2 != e; it2++){
                        cout << "" << *it2 << " ";
                        for (auto it3 = next(it2); it3 != e; it3++){
                            pair<int, int> p;
                            p.first = *it2;
                            p.second = *it3;
                            candidatos[p] = 0;
                            cont_pares_candidatos++;
                        }
                    }
                }
            }
            cout << "\n" << "Numero de pares candidatos en banda 2: " << cont_pares_candidatos << endl;
            cont_pares_candidatos = 0;
        }
    }
    
    cout << "Cantidad de candidatos: " << candidatos.size() << endl;
    std::map<pair<int, int>, double>::iterator it = candidatos.begin();
    cout << "Calculo de similitud de Jaccard..." << endl;
    while (it != candidatos.end()){
        it->second = jaccardSimilarity(it->first.first, it->first.second);
        it++;
    }
} 

void calcularSimilitud(){
    cout << "VALOR DE THRESHOLD: " << t << endl;
    int similares = 0;
    std::map<pair<int, int>, double>::iterator it = candidatos.begin();
    while (it != candidatos.end()){
        /* cout << " " << it -> second; */
        if (it->second >= t){
            cout << "Similitud(" << (it->first.first) + 1 << ", " << (it->first.second) + 1 << ") = " << it->second << endl;
            similares++;
        }
        it++;
    }
    if (similares == 0){
        cout << "No hay proteinas significativamente similares sobre ese umbral." << endl;
    }
}

int main(){
    ifstream inFile;
    inFile.open("proteinasunicas.fasta");
    string line;
    int line_length;
    bool flag=false;
    auto start = high_resolution_clock::now();
    srand(time(NULL));
    cout << "Leyendo archivo..." << endl;
    while(inFile){ 
        getline(inFile, line);
        line_length=line.length();
        /* cout << (doc_cont + 1) << endl; */
        doc_cont++;
        if (line.size() == 0){
            break;
        }
        for (int i = 0; i <= line_length-k; i++){
            string temp_string = line.substr(i, k);
            /* cout << temp_string << endl; */
            k_shingles.push_back(temp_string);
        }
        /////  i máximo define la dimension de cada signature representativa.
        for (int i = 0; i < signature_size; i++){
            /// Minhashing
            minhashing = Shingles().compute(k_shingles);
            vector_signature.push_back(minhashing);
        }
        /////  Signature matrix para entregar a LSH
        signature_matrix.push_back(vector_signature);
        k_shingles.clear();
        vector_signature.clear();
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Archivo leido y creada la matriz de signatures.-" << endl;
    cout << "Lineas procesadas: " << doc_cont << endl;
    cout << "Signature matrix size: " << signature_matrix.size() << " columnas.-" << endl;
    cout << "Duración de lectura de archivo y construccion de la matriz de signatures: " << (duration.count())/60 << " minutos, "<< (duration.count())%60 << " segundos.\n";
    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    auto start2 = high_resolution_clock::now();
    generarCandidatos();
    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<seconds>(stop2 - start2);
    cout << "Bandas: "<< bandas << endl;
    cout << "Filas por Banda: "<< filas_por_banda << endl;
    cout << "Duración de la generacion de candidatos para LSH: " << (duration2.count())/60 << " minutos, "<< (duration2.count())%60 << " segundos.\n";
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    auto start3 = high_resolution_clock::now();
    calcularSimilitud();
    auto stop3 = high_resolution_clock::now();
    auto duration3 = duration_cast<seconds>(stop3 - start3);
    cout << "Duración del calculo de la similitud entre los candidatos: " << (duration3.count())/60 << " minutos, "<< (duration3.count())%60 << " segundos.\n";
    ////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////
    inFile.close();
    return 0;
}
