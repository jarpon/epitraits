/*!
 * \class  CDFTools
 * \author Philippe Andrey, INRA
 * \date   2009.11.10 - création (pa)
 * \brief  Manipulation de fonctions de répartition
****************************************************************/

#include "cdftools2.h"

#include <usageerror.h>

#include <cmath>

 #define TRACE
#include "trace.h"

template<typename T>
CDFTools<T>::CDFTools()
{
}

template<typename T>
CDFTools<T>::~CDFTools()
{
}

/*! Calcule la fonction de répartition de l'échantillon de valeurs stockées
 * dans \c x.
 *
 * Le vecteur retourné donne les valeurs de la fonction de
 * répartition aux positions stockées dans \c x.
 *
 * \pre Le vecteur \c x doit être trié en ordre croissant.
 *
 * \exception UsageError Le vecteur \c x n'est pas trié en ordre croissant.
****************************************************************/
template<typename T>
Vector<T> CDFTools<T>::cdf(const Vector<T>& x) const
{
  if ( x.isIncreasing() == false )
  {
    UsageError usageError;
    usageError.setWhere( "Vector<T> CDFTools<T>::cdf(const Vector<T>&) const" );
    usageError.setWhat( "The input vector is not sorted in increasing order" );
    throw usageError;
  }

  const int n = x.getSize();
  Vector<T> y( n );

  for (int i = 0; i < n; ++i)
  {
    y[i] = T(i+1) / n;
  }

  return y;
}

/*! Calcule la fonction de répartition de l'échantillon de valeurs stockées
 * dans \c x.
 *
 * Le vecteur retourné donne les valeurs de la fonction de
 * répartition aux positions données par \c xEvals.
 *
 * \pre Le vecteur \c x doit être trié en ordre croissant.
 *
 * \exception UsageError Le vecteur \c x n'est pas trié en ordre croissant.
****************************************************************/
template<typename T>
Vector<T> CDFTools<T>::cdf(const Vector<T>& x, const Vector<T>& xEvals) const
{
  if ( x.isIncreasing() == false )
  {
    UsageError usageError;
    usageError.setWhere( "Vector<T> CDFTools<T>::cdf(const Vector<T>&, const Vector<T>&) const" );
    usageError.setWhat( "The input vector is not sorted in increasing order" );
    throw usageError;
  }

  const int n = x.getSize();
  const int numEvals = xEvals.getSize();
  Vector<T> y( numEvals );
  int i, p;

  for (i = 0, p = 0; i < n; ++i)
  {
    while ( p < numEvals && xEvals[p] < x[i] )
    {
      y[p++] = T(i) / n;
    }
  }

  while ( p < numEvals )
  {
    y[p++] = 1.0;
  }

  return y;
}

/*! Retourne la moyenne des fonctions de répartition d'un ensemble
 * d'échantillons de valeurs (non nécessairement tous de même taille).
 *
 * Pour chacun des échantillons stockés dans \c x, la méthode calcule
 * les valeurs de sa fonction de répartition aux positions données par
 * \c xEvals.
 * Retourne la moyenne de toutes les fonctions de répartition ainsi
 * calculées.
****************************************************************/
template<typename T>
Vector<T> CDFTools<T>::average(const vector<const Vector<T>*>& x, const Vector<T>& xEvals) const
{
  const int n = x.size();
  const int numEvals = xEvals.getSize();
  Vector<T> y( numEvals );

  y.setZeros();
  for (int i = 0; i < n; ++i)
  {
    y += cdf( *x[i], xEvals );
  }
  y /= n;

  return y;
}

/*! Retourne la moyenne des fonctions de répartition d'un ensemble
 * d'échantillons de valeurs (tous de même taille).
 *
 * Les échantillons sont donnés par les lignes de la matrice \c x.
 * La méthode calcule les valeurs de la fonction de répartition de chaque
 * échantillon aux positions données par \c xEvals.
 * Elle retourne la moyenne de toutes les fonctions de répartition ainsi
 * calculées.
****************************************************************/
template<typename T>
Vector<T> CDFTools<T>::average(const Matrix<T>& x, const Vector<T>& xEvals) const
{
  const int n = x.getSize1();
  vector<const Vector<T>*> v( n );

  for (int i = 0; i < n; ++i)
  {
    v[i] = &x[i];
  }

  return average( v, xEvals );
}

/*! Retourne la moyenne des fonctions de répartition d'un ensemble
 * d'échantillons de valeurs (non nécessairement tous de même taille).
 *
 * La méthode calcule la fonction de répartition de chacun des
 * échantillons stockés dans le vecteur \c x.
 * Les fonctions de répartitions sont évaluées sur l'ensemble des positions
 * données par \c x.
 * La méthode retourne la moyenne de toutes les fonctions de répartition
 * ainsi calculées.
****************************************************************/
template<typename T>
Vector<T> CDFTools<T>::average(const vector<const Vector<T>*>& x) const
{
  const int n = x.size();
  Vector<T> xEvals;

  for (int i = 0; i < n; ++i)
  {
    xEvals.append( *x[i] );
  }
  xEvals.sort();

  return average( x, xEvals );
}

/*! Retourne la moyenne des fonctions de répartition d'un ensemble
 * d'échantillons de valeurs (tous de même taille).
 *
 * Les échantillons sont donnés par les lignes de la matrice \c x.
 * La méthode calcule la fonction de répartition de chacun des
 * échantillons.
 * Les fonctions de répartitions sont évaluées sur l'ensemble des positions
 * données par \c x.
 * La méthode retourne la moyenne de toutes les fonctions de répartition
 * ainsi calculées.
****************************************************************/
template<typename T>
Vector<T> CDFTools<T>::average(const Matrix<T>& x) const
{
  Vector<T> xEvals;
  xEvals = x.cat();
  xEvals.sort();

  return average( x, xEvals );
}

/*! Pour convenance.
 *
 * Retourne l'enveloppe d'un percentile calculé sur un ensemble
 * d'échantillons de valeurs (non nécessairement tous de même taille).
 *
 * Les échantillons sont donnés par le vecteur \c x.
 * L'enveloppe est calculée aux positions données par les échantillons.
****************************************************************/
template<typename T>
Vector<T> CDFTools<T>::percentile(const vector<const Vector<T>*>& x, const float percent) const
{
  const int n = x.size();
  Vector<T> xEvals;

  for (int i = 0; i < n; ++i)
  {
    xEvals.append( *x[i] );
  }
  xEvals.sort();

  return percentile( x, xEvals, percent );
}

/*! Pour convenance.
 *
 * Retourne l'enveloppe d'un percentile calculé sur un ensemble
 * d'échantillons de valeurs (tous de même taille).
 *
 * Les échantillons sont donnés par les lignes de la matrice \c x.
 * L'enveloppe est calculée aux positions données par les échantillons.
****************************************************************/
template<typename T>
Vector<T> CDFTools<T>::percentile(const Matrix<T>& x, const float percent) const
{
  Vector<T> xEvals;
  xEvals = x.cat();
  xEvals.sort();

  return percentile( x, xEvals, percent );
}

/*! Pour convenance.
 *
 * Retourne l'enveloppe d'un percentile calculé sur un ensemble
 * d'échantillons de valeurs (non nécessairement tous de même taille).
 *
 * L'enveloppe est calculée aux positions données par \c xEvals.
****************************************************************/
template<typename T>
Vector<T> CDFTools<T>::percentile(
  const vector<const Vector<T>*>& x,
  const Vector<T>& xEvals,
  const float percent) const
{
  Matrix<T> matrix;
  Vectorf percents( 1 );
  percents[0] = percent;
  matrix = percentiles( x, xEvals, percents );
  return matrix[0];
/*  const int numSamples = x.size();
  const int numEvals = xEvals.getSize();
  Matrix<T> matrix( numEvals, numSamples );
  Vector<T> y( numEvals );

  for (int i = 0; i < numSamples; i++)
  {
    matrix.setColumn( i, cdf(*x[i],xEvals) );
  }

  for (int j = 0; j < numEvals; j++)
  {
    matrix[j].sort();
    y[j] = matrix[j].percentile( percent );
  }

  return y;*/
}

/*! Pour convenance.
 *
 * Retourne l'enveloppe d'un percentile calculé sur un ensemble
 * d'échantillons de valeurs (tous de même taille).
 *
 * Les échantillons sont donnés par les lignes de la matrice \c x.
 * L'enveloppe est calculée aux positions données par \c xEvals.
****************************************************************/
template<typename T>
Vector<T> CDFTools<T>::percentile(
  const Matrix<T>& x,
  const Vector<T>& xEvals,
  const float percent) const
{
  const int n = x.getSize1();
  vector<const Vector<T>*> v( n );

  for (int i = 0; i < n; ++i)
  {
    v[i] = &x[i];
  }

  return percentile( v, xEvals, percent );
}

/*! Retourne les enveloppes de percentiles calculés sur un ensemble
 * d'échantillons de valeurs (non nécessairement tous de même taille).
 *
 * Les échantillons sont donnés par le vecteur \c x.
 * Leurs fonctions de répartitions sont évaluées aux positions données
 * par \c xEvals.
 * En chaque position d'évaluation, la méthode détermine les valeurs
 * correspondant aux percentiles donnés par \c percents.
 * Les enveloppes ainsi calculées sont stockées ligne à ligne dans la
 * matrice renvoyée par la méthode.
 *
 * \pre Les valeurs données par \c percents doivent être comprises
 * entre 0 et 100.
 * La validité de cette condition n'est pas vérifiée par la méthode.
****************************************************************/
template<typename T>
Matrix<T> CDFTools<T>::percentiles(
  const vector<const Vector<T>*>& x,
  const Vector<T>& xEvals,
  const Vectorf& percents) const
{
  const int numSamples = x.size();
  const int numEvals = xEvals.getSize();
  const int numPercents = percents.getSize();
  Matrix<T> matrix( numEvals, numSamples );
  Matrix<T> y( numPercents, numEvals );

  for (int i = 0; i < numSamples; ++i)
  {
    matrix.setColumn( i, cdf(*x[i],xEvals) );
  }

  for (int j = 0; j < numEvals; ++j)
  {
    matrix[j].sort();
    for (int p = 0; p < numPercents; ++p)
    {
      y(p,j) = matrix[j].percentile( percents[p] );
    }
  }

  return y;
}

/*! Pour convenance.
 *
 * Retourne les enveloppes de percentiles calculés sur un ensemble
 * d'échantillons de valeurs (tous de même taille).
 *
 * Les échantillons sont donnés par les lignes de la matrice \c x.
 * Les enveloppes sont calculées aux positions données par \c xEvals.
****************************************************************/
template<typename T>
Matrix<T> CDFTools<T>::percentiles(
  const Matrix<T>& x,
  const Vector<T>& xEvals,
  const Vectorf& percents) const
{
  const int n = x.getSize1();
  vector<const Vector<T>*> v( n );

  for (int i = 0; i < n; ++i)
  {
    v[i] = &x[i];
  }

  return percentiles( v, xEvals, percents );
}

/*! Calcule les écarts entre les deux fonctions de répartition
 * empiriques \c cdf1 et \c cdf2.
 *
 * Les deux vecteurs donnent les points où des sauts ont lieu
 * dans les fonctions de répartition.
 * Les valeurs sont supposées avoir été triées au préalable en
 * ordre croissant.
 *
 * Trois écarts sont calculés:
 * - l'écart maximum, qu'il soit positif ou négatif
 * - l'écart maximum positif (càd max des Fcdf1 - Fcdf2)
 * - l'écart maximum négatif (càd min des Fcdf1 - Fcdf2)
 *
 * Chacun de ces écarts est retourné dans la seconde composante
 * du vecteur correspondant parmi \c maxDiff, \c maxDiffAbove et
 * \c maxDiffBelow.
 * La première composante de chacun de ces vecteurs donne la
 * position où l'écart maximum a été observé.
****************************************************************/
template<typename T>
void CDFTools<T>::differences(
  const Vector<T>& cdf1,
  const Vector<T>& cdf2,
  Vector<T>& maxDiff,
  Vector<T>& maxDiffAbove,
  Vector<T>& maxDiffBelow) const
{
  if ( cdf1.isIncreasing() == false || cdf2.isIncreasing() == false )
  {
    UsageError usageError;
    usageError.setWhere( "void CDFTools<T>::differences(...) const" );
    usageError.setWhat( "Non-increasing number sequence passed as cdf" );
    throw usageError;
  }

  const int n = cdf1.getSize();
  const int m = cdf2.getSize();
  int i = 0, j = 0, iprev, jprev;
  T diff;

  maxDiff.setZeros( 2 );
  maxDiffAbove.setZeros( 2 );
  maxDiffBelow.setZeros( 2 );

  while ( i < n && j < m )
  {
    iprev = i;
    jprev = j;

    if ( cdf1[i] <= cdf2[j] )
    {
      while ( i < n-1 && cdf1[i+1] == cdf1[i] ) i++;
      i++;
    }

    if ( cdf2[j] <= cdf1[iprev] )
    {
      while ( j < m-1 && cdf2[j+1] == cdf2[j] ) j++;
      j++;
    }

    diff = T(i)/n - T(j)/m;

    if ( diff > 0.0 && diff > maxDiffAbove[1] )
    {
      maxDiffAbove[0] = 0.5 * (cdf1[iprev]+cdf2[jprev]);
      maxDiffAbove[1] = diff;
    }

    if ( diff < 0.0 && diff < maxDiffBelow[1] )
    {
      maxDiffBelow[0] = 0.5 * (cdf1[iprev]+cdf2[jprev]);
      maxDiffBelow[1] = diff;
    }
  }

  maxDiff = maxDiffAbove[1] >= fabs(maxDiffBelow[1])? maxDiffAbove: maxDiffBelow;
}

template<typename T>
void CDFTools<T>::differences(
  const Vector<T>& cdf1,
  const Vector<T>& x1,
  const Vector<T>& cdf2,
  const Vector<T>& x2,
  Vector<T>& maxDiff,
  Vector<T>& maxDiffAbove,
  Vector<T>& maxDiffBelow) const
{
  if ( cdf1.isIncreasing() == false || cdf2.isIncreasing() == false )
  {
    UsageError usageError;
    usageError.setWhere( "void CDFTools<T>::differences(...) const" );
    usageError.setWhat( "Non-increasing number sequence passed as cdf" );
    throw usageError;
  }

  const int n = cdf1.getSize();
  const int m = cdf2.getSize();
  int i = 0, j = 0, iprev, jprev;
  T diff, v1 = 0, v2 = 0;

  maxDiff.setZeros( 2 );
  maxDiffAbove.setZeros( 2 );
  maxDiffBelow.setZeros( 2 );

  while ( i < n && j < m )
  {
    iprev = i;
    jprev = j;

    if ( cdf1[i] < cdf2[j] )
    {
      v1 = x1[i++];
    }
    else if ( cdf2[j] < cdf1[i] )
    {
      v2 = x2[j++];
    }
    else // les deux positions sont identiques
    {
      v1 = x1[i++];
      v2 = x2[j++];
    }

    diff = v1 - v2;

    if ( diff > 0.0 && diff > maxDiffAbove[1] )
    {
      maxDiffAbove[0] = 0.5 * (cdf1[iprev]+cdf2[jprev]);
      maxDiffAbove[1] = diff;
    }

    if ( diff < 0.0 && diff < maxDiffBelow[1] )
    {
      maxDiffBelow[0] = 0.5 * (cdf1[iprev]+cdf2[jprev]);
      maxDiffBelow[1] = diff;
    }
  }

  maxDiff = maxDiffAbove[1] >= fabs(maxDiffBelow[1])? maxDiffAbove: maxDiffBelow;
}


template<typename T>
Vector<T> CDFTools<T>::areasDifference(
  const Vector<T>& y1,
  const Vector<T>& x1,
  const Vector<T>& y2,
  const Vector<T>& x2) const
{
  if ( y1.isIncreasing() == false || y2.isIncreasing() == false )
  {
    UsageError usageError;
    usageError.setWhere( "void CDFTools<T>::areasDifference(...) const" );
    usageError.setWhat( "The cdf passed is not an increasing number sequence" );
    throw usageError;
  }

  if ( x1.getSize() != y1.getSize() || x2.getSize() != y2.getSize() )
  {
    UsageError usageError;
    usageError.setWhere( "void CDFTools<T>::areasDifference(...) const" );
    usageError.setWhat( "The sizes of at least one cdf passed and its repartition function do not match to each other" );
    usageError.setWhat( "Use Vector<T> CDFTools<T>::cdf(const Vector<T>&) if needed..." );
    throw usageError;
  }

  Vector<T> areasDifference;

  //   |   __|
  //   |   |   -->y[0] -- y axis
  //   |  _|
  //   |_|____ -->x[0]   -- x axis
  //

//  EVAL(y1);
//  EVAL(x1);
//  EVAL(y2);
//  EVAL(x2);

  T area1, area2;
  area1 = 0;
  area2 = 0;

  for ( int i = 1; i < y1.getSize(); ++i )
    area1 += ( x1[i]-x1[i-1] ) * y1[i-1];

  EVAL(area1);

  for ( int i = 1; i < y2.getSize(); ++i )
    area2 += ( x2[i]-x2[i-1] ) * y2[i-1];

  EVAL(area2);

  // output with 4 indexes
  areasDifference.setSize( 4 );

  // 1) Difference of areas
  areasDifference[0] = abs( area1 - area2 );
  // 2) Difference of areas ^ 2
  areasDifference[1] = abs( pow( area1, 2.) - pow( area2, 2.) );
  // 3) Coefficient between areas
  if ( area1 < area2 ) areasDifference[2] = area1 / area2;
  else areasDifference[2] = area2 / area1;
  // 4) Coefficient between areas ^ 2
  if ( area1 < area2 ) areasDifference[3] = pow( area1, 2.) / pow( area2, 2.);
  else areasDifference[3] = pow( area2, 2.) / pow( area1, 2.);

  return areasDifference;
}


template<typename T>
void CDFTools<T>::differences(const Vector<T>& cdf1, const Vector<T>& cdf2, Vector<T>& maxDiff) const
{
  Vector<T> maxDiffAbove;
  Vector<T> maxDiffBelow;

  differences( cdf1, cdf2, maxDiff, maxDiffAbove, maxDiffBelow );
}

/*! Positionne une fonction de répartition parmi un ensemble de fonctions
 * de répartition.
 *
 * La méthode calcule la fonction de répartition de l'échantillon de valeurs
 * \c x ainsi que celle de chacun des échantillons stockés ligne à ligne dans
 * la matrice \c distributionSamples.
 * Elle calcule ensuite l'écart entre chaque fonction de répartion et la
 * fonction de répartition moyenne \c distributionAverage (non nécessairement
 * calculée sur \c distributionSamples).
 * Elle retourne le rang de l'écart de l'échantillon \c x parmi tous les
 * écarts ainsi calculés.
****************************************************************/
template<typename T>
int CDFTools<T>::rank(
  const Vector<T>& x,
  const Matrix<T>& distributionSamples,
  const Vector<T>& xEvals,
  const Vector<T>& distributionAverage) const
{
  const int numSamples = distributionSamples.getSize1();
  Vector<T> maxDifferences( 1+numSamples );
  Vector<T> maxDiffBuffer;
  Vector<T> ignored1, ignored2;
  T xMaxDiff;

  differences( x, cdf(x), xEvals, distributionAverage, maxDiffBuffer, ignored1, ignored2 );
  xMaxDiff = maxDifferences[0] = maxDiffBuffer[1];
  for (int i = 0; i < numSamples; ++i)
  {
    differences(
      distributionSamples[i], cdf(distributionSamples[i]),
      xEvals, distributionAverage,
      maxDiffBuffer, ignored1, ignored2 );
    maxDifferences[i+1] = maxDiffBuffer[1];
  }
  maxDifferences.sort();

  int r = 0;
  while ( maxDifferences[r] != xMaxDiff )
  {
    ++r;
  }

  return r;
}

/*! Positionne une fonction de répartition parmi un ensemble de fonctions
 * de répartition.
 *
 * La méthode calcule la fonction de répartition de l'échantillon de valeurs
 * \c x ainsi que celle de chacun des échantillons stockés ligne à ligne dans
 * la matrice \c distributionSamples.
 * Elle calcule ensuite l'écart entre chaque fonction de répartion et la
 * fonction de répartition moyenne \c distributionAverage (non nécessairement
 * calculée sur \c distributionSamples).
 * Elle retourne le rang de l'écart de l'échantillon \c x parmi tous les
 * écarts ainsi calculés.
****************************************************************/
template<typename T>
Vector<float> CDFTools<T>::rankAndMaxDiff(
  const Vector<T>& x,
  const Matrix<T>& distributionSamples,
  const Vector<T>& xEvals,
  const Vector<T>& distributionAverage) const
{
  const int numSamples = distributionSamples.getSize1();
  Vector<T> maxDifferences( 1+numSamples );
  Vector<T> maxDiffBuffer;
  Vector<T> ignored1, ignored2;
  T xMaxDiff;

  differences( x, cdf(x), xEvals, distributionAverage, maxDiffBuffer, ignored1, ignored2 );
  xMaxDiff = maxDifferences[0] = maxDiffBuffer[1];
  for (int i = 0; i < numSamples; ++i)
  {
    differences(
      distributionSamples[i], cdf(distributionSamples[i]),
      xEvals, distributionAverage,
      maxDiffBuffer, ignored1, ignored2 );
    maxDifferences[i+1] = maxDiffBuffer[1];
  }
  maxDifferences.sort();

  int r = 0;
  while ( maxDifferences[r] != xMaxDiff )
  {
    ++r;
  }

  Vector<float> output( 2 );
  output[0] = r;
  output[1] = xMaxDiff;

  return output;
}

template class CDFTools<float>;
template class CDFTools<double>;
template class CDFTools<long double>;
