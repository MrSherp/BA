#ifndef IMAGELOADSAVE_H
#define IMAGELOADSAVE_H

#include <cstdio>
#include <bitmap_image.hpp>

#include "core.h"

template < typename ImageType >
void loadColorBitmap ( const std::string Filename, ImageType& ImageRed, ImageType& ImageGreen, ImageType& ImageBlue ) {

  bitmap_image rawImage ( Filename.c_str () );

  if ( !rawImage )
    throw std::runtime_error ( "Error while loading image!" );

  unsigned char data[3];

  const unsigned int height = rawImage.height ();
  const unsigned int width = rawImage.width ();

  ImageRed.resize ( width, height );
  ImageGreen.resize ( width, height );
  ImageBlue.resize ( width, height );

  for ( std::size_t y = 0; y < height; ++y ) {
    for ( std::size_t x = 0; x < width; ++x ) {
      rawImage.get_pixel ( x, y, data[0], data[1], data[2] );
      ImageRed ( x, y ) = static_cast < int > ( data[0] );
      ImageGreen ( x, y ) = static_cast < int > ( data[1] );
      ImageBlue ( x, y ) = static_cast < int > ( data[2] );
    }
  }
}

template < typename ImageType >
void saveColorBitmap ( const std::string Filename, const ImageType& ImageRed, const ImageType& ImageGreen, const ImageType& ImageBlue ) {

  const unsigned int width = ImageRed.rows ();
  const unsigned int height = ImageRed.cols ();

  bitmap_image rawImage ( width, height );

  for ( std::size_t y = 0; y < height; ++y ) {
    for ( std::size_t x = 0; x < width; ++x ) {
      rawImage.set_pixel ( x, y, ImageRed ( x, y ), ImageGreen ( x, y ), ImageBlue ( x, y ) );
    }
  }
  if ( boost::filesystem::exists ( Filename ) ) {
    co::Core::displayWarning ( "File already exists!" );
    std::remove ( Filename.c_str () );
  }
  rawImage.save_image ( Filename.c_str () );
}

template < typename ImageType >
void saveBitmap ( const std::string Filename, const ImageType& Image ) {
  saveColorBitmap ( Filename, Image, Image, Image );
}

template < typename ImageType >
void convertToGrayScale ( const ImageType& ImageRed, const ImageType& ImageGreen, const ImageType& ImageBlue, ImageType& Dest ) {
  const unsigned int xMax = ImageRed.rows ();
  const unsigned int yMax = ImageRed.cols ();

  Dest.resize ( xMax, yMax );

  for ( std::size_t y = 0; y < yMax; ++y ) {
    for ( std::size_t x = 0; x < xMax; ++x ) {
      Dest ( x, y ) = static_cast < int > ( sqrt ( co::sqr ( ImageRed ( x, y ) ) + co::sqr ( ImageGreen ( x, y ) ) + co::sqr ( ImageBlue ( x, y ) ) ) / static_cast < long double > ( 3 ) );
    }
  }
}

int addCounterToSaveDirectory ( const std::string CounterFileName ) {
  if ( !boost::filesystem::exists( CounterFileName ) ) {
    std::ofstream out ( CounterFileName.c_str() );
    out << 0 << std::endl;
    out.close ( );
  }
  std::fstream counterFile;
  counterFile.open ( CounterFileName.c_str () );
  int counter = 0;
  if ( counterFile.is_open () ) {
    std::string temp;
    std::getline ( counterFile, temp );
    counter = atoi ( temp.c_str () );
    counterFile.seekg ( std::ios::beg );
    counterFile << ++counter;
  }
  else
    throw std::runtime_error ( "Cannot open parameter file for writing!" );
  counterFile.close ();
  return counter;
}

#endif
