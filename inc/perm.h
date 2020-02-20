#include <iostream>
#include <stdlib.h>

using namespace std;

void printlist(int* list, int size)
{
  for (int i=0;i<size;i++) cout << list[i];
  cout <<endl;
}
void swap(int *list,int i1, int i2) // swap two elements of list
{
  int tmp = list[i2];
  list[i2] = list[i1];
  list[i1] = tmp;
}

void sort(int * list, int size) // numerically sort the list
{
  for (int i=0;i<size;i++)
    if (list[i] < list[0]) swap(list,0,i);
  if (size>1) sort(&list[1],size-1);
}
void get_nth_perm(int *list, int size, int n)
{
  if (size <=1 || n==0) return;
  sort(list,size);

  int i,X=1;
  for(i=1;i<size;i++) X*=i;

  for (i=n/X;i>0;i--) swap(list, i,i-1);

  get_nth_perm(&list[1],size-1,n%X);
}
bool get_next_perm(int *list, int size)
{
  int i,j;
  if (size < 2) return false; // false means done
  for(i=size-1; i>0 && list[i] < list[i-1]; i--);
  if (i==0) return false; // false means done (already on last element)
  sort(&list[i],size-i);
  i--;
  for(j=i+1;j<size && list[i]>list[j]; j++);
  swap(list,i,j);
  return true; // true means there are more permutations to be had
}
