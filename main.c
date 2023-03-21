#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "string.h"

#define M 4
#define OPTIONS 8
#define VISITED -1
#define INVALID -1
#define N 6
#define LETTERS 26
#define EMPTY 0


bool check (int hist1[N], int hist2[N]);


bool check (int hist1[N], int hist2[N]){
    for (int i=0; i<N; i++){
        if ((hist1[i]>0 && hist2[i]==0) || (hist2[i]>0 && hist1[i]==0)){
            return false;
        }
    }
    return true;
}

bool hist_compare(int hist1[], int hist2 []){
    for(int i = 0; i<N; i++){
        if(hist1[i]!= hist2[i]){
            return false;
        }
    }
    return true;
}

int str_length(char* str){
    int len = 0;
    while(*str != '\0'){
        if(*str != ' ')
            len++;
        str++;
    }
    return len;
}
int count_mixed_str(char* s1, char* s2){
    int hist1[N] = {0};
    int hist2[N] = {0};
    int counter = 0;
    int len1 = str_length(s1);
    int len2 = str_length(s2);
    for(int i = 0; i<len1; i++){
        hist1[s1[i] - 'a']++;
        hist2[s2[i] - 'a']++;
    }
    for(int i = 0; i <= len2-len1; i++){
        if(hist_compare(hist1, hist2)){
            counter++;
        }
        hist2[s2[i]-'a']--;
        hist2[s2[len1+i] - 'a']++;
    }
    return counter;
}



int num (char* s);
bool legal (char *str, int n);

bool legal(char *str, int n){
    int counter_int = 0, counter_big_letter=0, counter_small_letter=0;
    for(int i=0; i<n; i++){
        if ((str[i]<='9') && (str[i]>='0')){
            counter_int++;
        }
        if ((str[i]<='Z') && (str[i]>='A')){
            counter_big_letter++;
        }
        if ((str[i]<='z') && (str[i]>='a')){
            counter_small_letter++;
        }
    }
    if (((counter_big_letter==1) && (counter_small_letter==4) && (counter_int==1)) || ((counter_big_letter<=3) && (counter_int==3))){
        return true;
    }
    else{
        return false;
    }
}

int num (char* s){
    int counter=0;
    int length = strlen (s);
    for(int i=0; i<=length-3; i++){
        for(int j=i+1; s[j]!='\0'; j++){
            if (legal (s,  j+1)==1){
                counter++;
            }
        }
        s++;
    }
    return counter;
}

int countChars(char * str, int n, char c){
    int low = 0, high = n-1, mid = 0;
    int first_index = 0, last_index = 0;
    while(low < high){
        mid = (low+high)/2;
        if(str[mid] == c){
            break;
        }
        if(str[mid] > c){
            high = mid-1;
        }
        else{
            low = mid +1;
        }
    }
    low = 0, high = mid;
    int temp_mid = 0;
    while(low <= high){
        temp_mid = (low+high)/2;
        if(str[temp_mid] == c && str[temp_mid] > str[temp_mid-1]){
            first_index = temp_mid;
            break;
        }
        if(str[temp_mid]!= c){
            low = temp_mid + 1;
        }
        else{
            high = temp_mid - 1;
        }
    }
    low = mid, high = n-1;
    while(low <= high){
        temp_mid = (low+high)/2;
        if((str[temp_mid] == c) && (str[temp_mid] < str[temp_mid+1]) || (low == high)){
            last_index = temp_mid;
            break;
        }
        if(str[temp_mid] == c){
            low = temp_mid + 1;
        }
        else{
            high = temp_mid - 1;
        }
    }
    return last_index - first_index + 1;
}

void CycleShift(int arr[], int n, int k){
    for(int i = 0; i<n-k; i+=k) {
        for(int j = 0, temp = 0; j<k; j++) {
            temp = arr[k+i+j];
            arr[k+i+j] = arr[j];
            arr[j] = temp;
        }
    }
}

void make_hist_from_index(char* str, int index, int m, int hist[]){
    for(int i = index; i < index + m; i++){
        hist[str[i]-'a']++;
    }
}

void clear_hist(int hist[]){
    for(int i = 0; i<26; i++) {
        hist[i] = 0;
    }
}

int CountPermutations(char* s1, char* s2){
    int hist1 [N] = {0}, hist2 [N] = {0};
    int len2 = str_length(s2);
    int len1 = str_length(s1);
    int counter =  0;
    make_hist_from_index(s2, 0, len2, hist2);
    for(int i = 0; i <= len1 - len2; i++){
        make_hist_from_index(s1, i, len2, hist1);
        if(hist_compare(hist1, hist2)){
            counter++;
        }
        clear_hist(hist1);
    }
    return counter;
}

int find_last_prev(int arr[], int n, int x)
{
    int left = 0, right = n - 1;
    while (left <= right)
    {
        int mid = left + (right - left) / 2;
        if (arr[mid] >= x)
            right = mid - 1;
        if (arr[mid] < x)
            left = mid + 1;
    }
    return left - 1;
}


int countOccurences(int arr[], int n, int x)
{
    int prev = find_last_prev(arr, n, x);
    int last = find_last_prev(arr, n, x + 1);
    return last - prev;
}

void make_hist(char* str, int n, int hist []){
    for(int i =0; i<n; i++){
        hist[str[i] - 'a']++;
    }
}

char toLower(char c) {
    if(c == ' ')
        return c;
    if (c >= 'a' && c <= 'z')
        return c;
    return c - 'A' + 'a';
}

bool areAnagrams(char *s1, char* s2){
    int hist1[N] = {0}, hist2 [N] = {0};
    int len1 = str_length(s1), len2 = str_length(s2);
    if(len2 != len1)
        return false;
    for(int i = 0; i<len1; i++){
        s1[i] = toLower(s1[i]);
        s2[i] = toLower(s2[i]);
    }
    make_hist(s1, len1, hist1);
    make_hist(s2, len2, hist2);
    return hist_compare(hist1, hist2);
}


int Countpermutations (char* s1, char* s2){ // problem here - fix
    int m = str_length(s2);
    int n = str_length(s1);
    int counter=0;
    int hist1[N]= {0};
    int hist2[N]= {0};

    //histogram for s2
    for (int i=0; i<m; i++){
        hist2[s2[i]-'a']++;
    }
    //histogram for every combination of s1
    for (int i=0; i<=n-m; i++){
        bool flag = false; // added for check
        clear_hist(hist1); // added for correctness
        for (int j=0; j<m; j++){
            hist1[s1[j+i]-'a']++;
        }
        for (int k=0; k<N; k++){
            if (hist1[k]!=hist2[k]){
                flag = true; // added for correctness
                break;
            }
        }
        if(!flag) // added for correctness
            counter++;
    }
    return counter;
}



int meeting_point(int a[], int na, int b[], int nb)
{
    if (na < 1 || nb < 1)
        return -1;
    if (na == 1 || nb == 1) {
        if (a[0] == b[0]) {
            return 0;
        }
        return -1;
    }
    int mid = 0;
    if(na <= nb){
        mid = na - 1;
    }
    else{
        mid = nb-1;
    }
    mid /= 2;
    if (a[mid] <= b[mid])
        return meeting_point(a, mid + 1, b, mid + 1);
    mid++;
    int result = meeting_point(a + mid, na - mid, b + mid, nb - mid);
    if(result==-1){
        return -1;
    }
    return mid+result;
}

int length_of_longest_stunning_substr(char *str)
{
    int hist[256] = { 0 }, max_unique = 0;
    char *left = str - 1, *right = left;
    while (*++right)
    {
        if (right - left - 1 > max_unique)
            max_unique = right - left - 1;
        if (++hist[*right] > 1)
        {
            do
            {
                --hist[*++left];
            } while (*left != *right);
        }
    }
    if (right - left - 1 > max_unique)
        max_unique = right - left - 1;
    return max_unique;
}

bool IsConstructibleAux(int numbers[], int n, int result, int currentSum)
{
    if (currentSum == result)
        return true;
    if (n == 0)
        return false;
    if (IsConstructibleAux(numbers + 1, n - 1, result, currentSum + numbers[0]))
        return true;
    if (IsConstructibleAux(numbers + 1, n - 1, result, currentSum - numbers[0]))
        return true;
    return IsConstructibleAux(numbers + 1, n - 1, result, currentSum);
}

bool IsConstructible(int numbers[],int n, int result)
{
    return IsConstructibleAux(numbers, n, result, 0);
}

bool islegal(int r, int c) { return (r >= 0 && r < N && c >= 0 && c < N); }

void disable(bool matrix[N][N], int r, int c) {
  if (!islegal(r, c) || matrix[r][c] == false)
    return;
  matrix[r][c] = false;
  disable(matrix, r, c - 1);
  disable(matrix, r, c + 1);
  disable(matrix, r - 1, c);
  disable(matrix, r + 1, c);
}

int cntTrueRegions(bool matrix[N][N]) {
  int counter = 0;
  for (int r = 0; r < N; r++) {
    for (int c = 0; c < N; c++) {
      if (matrix[r][c] == true) {
        counter++;
        disable(matrix, r, c);
      }
    }
  }
  return counter;
}

bool in_bounds(int n, int i, int j) {
    return i >= 0 && j >= 0 && i < n && j < M;
}

bool is_valid(char mat[][M], int i, int j, int next_i, int next_j,
              char letter) {
    if (0 > next_i || 0 > next_j || next_i > N || next_j > M) {
        return false;
    }
    if (mat[next_i][next_j] != letter) {
        return false;
    }
    return true;
}


int find_word_from_location(char mat[][M], int i, int j, char *word) {
    if (!(*word)) {
        return 1;
    }
    int count = 0;
    int rowNum[] = {-1, -1, -1, 0, 0, 1, 1, 1};
    int colNUm[] = {-1, 0, 1, -1, 1, -1, 0, 1};
    for (int dir = 0; dir < OPTIONS; dir++) {
        int next_i = i + rowNum[dir], next_j = j + colNUm[dir];
        if (is_valid(mat, i, j, next_i, next_j, word[0])) {
            mat[next_i][next_j] = VISITED;
            count += find_word_from_location(mat, next_i, next_j, word + 1);
            mat[next_i][next_j] = word[0];
        }
    }
    return count;
}

int find_num_of_occurences(char letters[N][M], char *word) {
    if ((*word) == '\0') {
        return 0;
    }
    int count = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            if (letters[i][j] == word[0]) {
                letters[i][j] = VISITED;
                count += find_word_from_location(letters, i, j, word + 1);
                letters[i][j] = word[0];
            }
        }
    }
    return count;
}

void findSubArray(int cumsum[], int n, int sum, int indices[2]){
    int first = 0, second=0, left_sum = 0;
    indices[0] = -1;
    indices[1] = -1;
    while(first < n){
        if(cumsum[first] - cumsum[second] > sum){
            second++;
        }
        else if(cumsum[first] - cumsum[second] < sum){
            first++;
        }
        else{
            indices[0] = second;
            indices[1] = first;
            return;
        }
    }
}
void swap(int*x, int *y){
    int temp = *x;
    *x = *y;
    *y = temp;
}
void stable_shuffle(int arr[], int n) {
    if (n <= 2) return;
    stable_shuffle(arr, n / 2);
    stable_shuffle(arr + n / 2, n / 2);
    for (int i = 0; i < n / 4; i++) {
        swap(arr + n / 4 + i, arr + n / 2 + i);
    }
}

int findPartiallySorted(int* arr, int n, int x) {
    int low = 0, high = n - 1;
    int mid, temp;
    while (low <= high) {
        mid = (low + high) / 2;
        if (arr[mid] == x)
            return mid;
        if (arr[mid] > x) {
            temp = mid;
            while (temp < high && arr[temp] > arr[temp + 1]) {
                temp++;
                if (arr[temp] == x)
                    return temp;
            }
            high = mid - 1;
        }
        if (arr[mid] < x) {
            temp = mid;
            while (temp > low && arr[temp] < arr[temp - 1]) {
                temp--;
                if (arr[temp] == x)
                    return temp;
            }
            low = mid + 1;
        }
    }
    return -1;
}


int findShortestMultiSub(char* s, char* t) {
    int len1 = strlen(s);
    int len2 = strlen(t);
    int hash_of_t[26] = {0};
    int hash_of_s[26] = {0};
    for (int i = 0; i < len2; i++)
        hash_of_t[t[i] - 'a']++;
    int start = 0, start_index = -1, min_len = len1;
    int count = 0;
    for (int j = 0; j < len1; j++) {
        hash_of_s[s[j] - 'a']++;
        if (hash_of_s[s[j] - 'a'] <= hash_of_t[s[j] - 'a'])
            count++;
        if (count == len2) {
            while (hash_of_s[s[start] - 'a'] > hash_of_t[s[start] - 'a'] ||
                   hash_of_t[s[start] - 'a'] == 0) {
                if (hash_of_s[s[start] - 'a'] > hash_of_t[s[start] - 'a'])
                    hash_of_s[s[start] - 'a']--;
                start++;
            }
            int len_window = j - start + 1;
            if (min_len > len_window) {
                min_len = len_window;
                start_index = start;
            }
        }
    }
    if (start_index == -1) {
        return -1;
    }
    return min_len;
}


int longest_path_aux(int roads[N][N], int src, int dest, bool *visited)
{
    int i = 0, max = -1;
    if(visited[src] == true)
        return INVALID;
    if(src == dest)
        return 1;
    visited[src] = true;
    for(i = 0; i < N; i++)
    {
        if(roads[src][i])
        {
            int k = longest_path_aux(roads, i, dest, visited);
            if(k+roads[src][i] > max){
                max = k+roads[src][i];
            }
        }
    }
    visited[src] = false;
    return max;
}

int get_max_letter(int hist[LETTERS], int last){
    int max_letter = -1;
    int max_counter = 0;
    for(int letter = 0; letter<LETTERS; letter++){
        if(hist[letter] > max_counter && letter != last){
            max_letter = letter;
            max_counter = hist[letter];
        }
    }
    return max_letter;
}

bool q3(char* s) {
    int hist[LETTERS] = {0};
    make_hist(s,  str_length(s), hist);
    int last = -1;
    for (int i = 0; s[i] != '\0'; i++) {
        int letter = get_max_letter(hist, last);
        if (letter == -1) {
            return false;
        }
        hist[letter]--;
        s[i] = letter + 'a';
        last = letter;
    }
    return true;
}

int longest_path(int roads[N][N], int src, int dest)
{
    bool visited[N] = {0};
    return longest_path_aux(roads, src, dest, visited);
}

void BuildHistogram(char* str, int* histogram)
{
    while (*str)
    {
        ++histogram[*str - 'A'];
        ++str;
    }
}
void permute(int letters[LETTERS], char* str, int length) {
    str[length] = '\0';
    printf("%s\n", str);
    for (int letter = 0; letter < LETTERS; ++letter) {
        if (letters[letter]) {
            --letters[letter];
            str[length] = 'A' + letter;
            permute(letters, str, length + 1);
            ++letters[letter];
        }
    }
}

void PrintAllPossibleSubstrings(char *str)
{
    int letters[LETTERS] = { 0 };
    char* word = malloc((strlen(str) + 1) * sizeof(char));
    BuildHistogram(str, letters);
    permute(letters, word, 0);
    free(word);
}

bool AreAllSumsEqual(int sums[], int k)
{
    for (int i = 1; i < k; ++i)
        if (sums[i] != sums[i - 1])
            return false;
    return true;
}

bool IsKSplittable_aux(int arr[], int n, int sums[], int k, int current)
{
    if (current == n) /// if all numbers used check if sums are right
        return AreAllSumsEqual(sums, k);
    for (int sum = 0; sum < k; sum++) /// we try to add to each k-sum
    {
        sums[sum] += arr[current]; // add it to the set
        if (IsKSplittable_aux(arr, n, sums, k, current + 1))
            return true;
        sums[sum] -= arr[current];
    }
    return false;
}

bool IsKSplittable(int arr[], int n, int k)
{
    int* sums = malloc(k * sizeof(int)); /// init the sum array
    if(!sums)
        return false;
    for (int i = 0; i < k; ++i)
        sums[i] = 0;
    bool isPossible = IsKSplittable_aux(arr, n, sums, k, 0);
    free(sums);
    return isPossible;
}


int main (){
    PrintAllPossibleSubstrings("");
    return 0;
}


