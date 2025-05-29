#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map> 
#include <array>

#include "contigs.hpp"

using namespace std;

typedef int Terminus;
const Terminus START = 0;
const Terminus RCSTART = 1;
const Terminus END = 2;
const Terminus RCEND = 3;

class Overlap {
public:
    // Конструктор класса Overlap
    Overlap(ContigIndex contig_i, Terminus terminus_i,
            ContigIndex contig_j, Terminus terminus_j,
            int ovl_len) :
            contig_i(contig_i), terminus_i(terminus_i),
            contig_j(contig_j), terminus_j(terminus_j),
            ovl_len(ovl_len) {}

    // Поля класса
    ContigIndex contig_i; // индекс (ключ) первого контига
    Terminus terminus_i; // термин (первого контига) участвующий в перекрытии
    ContigIndex contig_j; // индекс (ключ) второго контига
    Terminus terminus_j; // термин (второго контига) участвующий в перекрытии
    int ovl_len; // длина перекрытия

    // Переопределение оператора преобразования в строку (аналог __repr__ в Python)
    std::string to_string() const {                                                 /////////Не используется (для тестов)
        return "<" + std::to_string(contig_i) + "-" + std::to_string(terminus_i) +
               "; " + std::to_string(contig_j) + "-" + std::to_string(terminus_j) +
               "; len=" + std::to_string(ovl_len) + ">";
    }

    // Переопределение оператора сравнения для проверки на равенство (аналог __eq__ в Python)
    bool operator==(const Overlap& other) const {
        return contig_i == other.contig_i &&
               terminus_i == other.terminus_i &&
               contig_j == other.contig_j &&
               terminus_j == other.terminus_j &&
               ovl_len == other.ovl_len;
    }

    // Переопределение оператора хеширования (аналог __hash__ в Python)
    size_t hash() const {
        size_t hash_value = 0;
        hash_combine(hash_value, contig_i);
        hash_combine(hash_value, terminus_i);
        hash_combine(hash_value, contig_j);
        hash_combine(hash_value, terminus_j);
        hash_combine(hash_value, ovl_len);
        return hash_value;
    }

private:
    // Вспомогательная функция для комбинирования хеша
    void hash_combine(size_t& seed, size_t value) const {
        seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
};

class OverlapCollection {
public:
    // Конструктор класса OverlapCollection
    OverlapCollection() {}

    // Метод для получения списка перекрытий, связанных с контигом по его ключу
    std::vector<Overlap> operator[](ContigIndex key) const {
        auto it = _collection.find(key);
        if (it != _collection.end()) {
            return it->second;
        } else {
            return std::vector<Overlap>(); // Возвращаем пустой вектор, если ключ не найден
        }
    }

    // Метод для получения размера коллекции
    size_t size() const {
        return _collection.size();
    }

    // Метод для добавления перекрытия в коллекцию
    void add_overlap(ContigIndex key, const Overlap& overlap) {
        _collection[key].push_back(overlap);
    }

    // Переопределение оператора преобразования в строку
    std::string to_string() const {                      /////////Не используется (для тестов)
        std::string result = "{";
        for (const auto& pair : _collection) {
            result += std::to_string(pair.first) + ": [";
            for (const auto& overlap : pair.second) {
                result += overlap.to_string() + ", ";
            }
            result += "], ";
        }
        result += "}";
        return result;
    }

        // Добавим методы begin и end для использования в цикле for
    auto begin() const { return std::begin(_collection); }
    auto end() const { return std::end(_collection); }

private:
    // Словарь, хранящий списки перекрытий для каждого контига
    std::unordered_map<ContigIndex, std::vector<Overlap>> _collection;
};

int find_overlap_s2s(const std::string& seq1, const std::string& seq2, int mink, int maxk) {
    maxk = std::min(maxk, std::min(static_cast<int>(seq1.length()), static_cast<int>(seq2.length())));
    if (maxk < mink) return 0;

    int overlap = 0;
    for (int i = mink; i <= maxk; ++i) {
        if (std::equal(seq1.begin(), seq1.begin() + i, seq2.begin())) {
            overlap = i;
        }
    }
    return overlap;
}

int find_overlap_e2s(const std::string& seq1, const std::string& seq2, int mink, int maxk) {
    maxk = std::min(maxk, std::min(static_cast<int>(seq1.length()), static_cast<int>(seq2.length())));
    if (maxk < mink) return 0;

    int overlap = 0;
    for (int i = mink; i <= maxk; ++i) {
        if (std::equal(seq1.end() - i, seq1.end(), seq2.begin(), seq2.begin() + i)) {
            overlap = i;
        }
    }
    return overlap;
}

int find_overlap_e2e(const std::string& seq1, const std::string& seq2, int mink, int maxk) {
    maxk = std::min(maxk, std::min(static_cast<int>(seq1.length()), static_cast<int>(seq2.length())));
    if (maxk < mink) return 0;

    int overlap = 0;
    for (int i = mink; i <= maxk; ++i) {
        if (std::equal(seq1.end() - i, seq1.end(), seq2.end() - i, seq2.end())) {
            overlap = i;
        }
    }
    return overlap;
}

OverlapCollection detect_adjacent_contigs(const ContigCollection& contig_collection,
                                          int mink, int maxk) {
    OverlapCollection overlap_collection;
    int num_contigs = contig_collection.size();
    
    #pragma omp parallel for schedule(dynamic) // Enable OpenMP parallelization
    for (ContigIndex i = 0; i < num_contigs; i++) {
        if (contig_collection[i].length <= mink) {
            std::cout << "\r" << i + 1 << "/" << num_contigs;
            continue;
        }

        std::vector<Overlap> local_overlaps;

        // Check self-overlaps first
        int ovl_len = find_overlap_e2s(contig_collection[i].end,
                                     contig_collection[i].start,
                                     mink, maxk);
        if (ovl_len > 0 && ovl_len < contig_collection[i].length) {
            local_overlaps.emplace_back(i, END, i, START, ovl_len);
            local_overlaps.emplace_back(i, START, i, END, ovl_len);
        }

        ovl_len = find_overlap_s2s(contig_collection[i].start,
                                  contig_collection[i].rcend,
                                  mink, maxk);
        if (ovl_len != 0) {
            local_overlaps.emplace_back(i, START, i, RCEND, ovl_len);
            local_overlaps.emplace_back(i, RCEND, i, START, ovl_len);
        }

        // Compare with other contigs
        for (ContigIndex j = i + 1; j < num_contigs; j++) {
            // Pre-calculate all possible overlaps for this pair
            std::array<int, 8> overlaps = {
                find_overlap_e2s(contig_collection[j].end, contig_collection[i].start, mink, maxk),
                find_overlap_e2s(contig_collection[i].end, contig_collection[j].start, mink, maxk),
                find_overlap_e2s(contig_collection[j].rcstart, contig_collection[i].start, mink, maxk),
                find_overlap_e2s(contig_collection[i].end, contig_collection[j].rcend, mink, maxk),
                find_overlap_s2s(contig_collection[i].start, contig_collection[j].start, mink, maxk),
                find_overlap_e2e(contig_collection[i].end, contig_collection[j].end, mink, maxk),
                find_overlap_s2s(contig_collection[i].start, contig_collection[j].rcend, mink, maxk),
                find_overlap_e2e(contig_collection[i].end, contig_collection[j].rcstart, mink, maxk)
            };

            // Add non-zero overlaps
            if (overlaps[0] != 0) {
                local_overlaps.emplace_back(i, START, j, END, overlaps[0]);
                local_overlaps.emplace_back(j, END, i, START, overlaps[0]);
            }
            if (overlaps[1] != 0) {
                local_overlaps.emplace_back(i, END, j, START, overlaps[1]);
                local_overlaps.emplace_back(j, START, i, END, overlaps[1]);
            }
            if (overlaps[2] != 0) {
                local_overlaps.emplace_back(i, START, j, RCSTART, overlaps[2]);
                local_overlaps.emplace_back(j, START, i, RCSTART, overlaps[2]);
            }
            if (overlaps[3] != 0) {
                local_overlaps.emplace_back(i, END, j, RCEND, overlaps[3]);
                local_overlaps.emplace_back(j, END, i, RCEND, overlaps[3]);
            }
            if (overlaps[4] != 0) {
                local_overlaps.emplace_back(i, START, j, START, overlaps[4]);
                local_overlaps.emplace_back(j, START, i, START, overlaps[4]);
            }
            if (overlaps[5] != 0) {
                local_overlaps.emplace_back(i, END, j, END, overlaps[5]);
                local_overlaps.emplace_back(j, END, i, END, overlaps[5]);
            }
            if (overlaps[6] != 0) {
                local_overlaps.emplace_back(i, START, j, RCEND, overlaps[6]);
                local_overlaps.emplace_back(j, RCEND, i, START, overlaps[6]);
            }
            if (overlaps[7] != 0) {
                local_overlaps.emplace_back(i, END, j, RCSTART, overlaps[7]);
                local_overlaps.emplace_back(j, RCSTART, i, END, overlaps[7]);
            }
        }

        // Critical section to add local overlaps to the shared collection
        #pragma omp critical
        {
            for (const auto& ovl : local_overlaps) {
                overlap_collection.add_overlap(ovl.contig_i, ovl);
            }
        }

        std::cout << "\r" << i + 1 << "/" << num_contigs;
    }

    std::cout << std::endl;
    return overlap_collection;
}