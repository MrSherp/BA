#ifndef __TIKZPLOTTER_H
#define __TIKZPLOTTER_H

//! run pdflatex for given tex-file
bool runPdflatex ( const string& CommandFileName, const string& SaveDirectory = "", const bool AddDotTex = false )  {
  string systemCommand = "pdflatex ";
  if ( SaveDirectory != "" ) systemCommand += " -output-directory " + SaveDirectory + " ";
  systemCommand += CommandFileName;
  if( AddDotTex ) systemCommand += ".tex";
  cout << systemCommand << endl;
  bool failed = true;
  for ( unsigned short int i = 0; i < 2; ++i ) {
    failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
    if ( failed ) cerr << "runPdflatex: Calling pdflatex returned an error.\n";
  }
  return !failed;
}




template< typename RealType >
class TikzPlotterHelperClass{
public :
    
    void generateIncludes ( std::ofstream &out,
         const RealType paperwidth, const RealType paperheight,
         const RealType leftBorder = 0.0, const RealType rightBorder = 0.0, const RealType topBorder = 0.0, const RealType bottomBorder = 0.0) const {

        out << "\\documentclass{article}" << endl;
        out << "\\usepackage[paperwidth=" << paperwidth << "cm,paperheight=" << paperheight << "cm]{geometry}" << endl;
        out << "\\geometry{"
            << "left=" << leftBorder << "cm," 
            << "right=" << rightBorder << "cm,"
            << "top=" << topBorder << "cm,"
            << "bottom=" << bottomBorder << "cm}" << endl
            << "\\usepackage[dvipsnames]{xcolor}" << endl
            << "\\usepackage{pdfpages}" << endl
            << "\\usepackage{pgfplots}" << endl
            << "\\pgfplotsset{compat=newest}" << endl
            << "\\pagestyle{empty}" << endl << endl
            << "\\usetikzlibrary{arrows}" << endl 
            << "\\usetikzlibrary{arrows.meta}" << endl << endl;
        // only for 3d plot
        out << "\\usepackage{tikz}" << endl
            << "\\usepackage{tikz-3dplot}" << endl
            << "\\usetikzlibrary{patterns}" << endl << endl << endl;
    }
    
    void generateBeginDocument ( std::ofstream &out ) const { out << "\\begin{document}" << endl << endl;}
    void generateEndDocument ( std::ofstream &out ) const { out << endl << endl << "\\end{document}" << endl;}
    void generateNewPage ( std::ofstream &out ) const { out << "\\newpage \\clearpage" << endl << endl; };
    
    void generateBeginTikzPicture ( std::ofstream &out, const RealType resizebox, const int dimension = 2, const string &comment = "" ) const {
        if ( comment != "" ) out << "% " << comment << endl;
        out << "%fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff" << endl
            << "\\begin{figure}[!htbp]" << endl
            << "\\resizebox{" << resizebox << "\\paperwidth}{!}{" << endl;
        switch( dimension ){
                case 1 :{
                    out << "\\begin{tikzpicture}[x=\\paperwidth,y=\\paperwidth]" << endl;
                }break;
            
                case 2 :{
                    out << "\\begin{tikzpicture}[x=\\paperwidth,y=\\paperwidth]" << endl;
                }break;
                
                case 3 :{
                    out << "\\tdplotsetmaincoords{70}{120} %rotation in x axis and z axis" << endl
                        << "\\begin{tikzpicture}[scale=10, tdplot_main_coords]" << endl;
                }break;
                
                default :
                        throw std::invalid_argument ( aol::strprintf ( "wrong mode in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
         }
    }
        
    void generateEndTikzPicture ( std::ofstream &out ) const {
        out << "\\end{tikzpicture}" << endl
            << "}" << endl
            << "\\end{figure}" << endl
            << "%fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff" << endl;
    }

    void generateBarChartAxis ( std::ofstream &out ) const {
        out << "\\draw (0,0) -- (1,0);" << endl
            << "\\draw (0,0) -- (0,-0.1);" << endl
            << "\\draw (1,0) -- (1,-0.1);" << endl
            << "\\draw (-0,0) -- (-0,1);" << endl << endl;   
    }
    
};








template< typename RealType, typename VectorType >
class TikzPlotterCurve{
private:
  
  TikzPlotterHelperClass<RealType> _tikzHelper;
  const string _saveDirectory;
  
public :
    
    TikzPlotterCurve ( const string saveDirectory ) : _saveDirectory ( saveDirectory ) { }
    
    void getMinMaxValues( const VectorType &vec, RealType &minVal, RealType &maxVal ) const {
        minVal = 2.e+19; maxVal = -2.e+19;
        for( int i=0; i<vec.size(); ++i ){
          if( vec[i] > maxVal ) maxVal = vec[i];
          if( vec[i] < minVal ) minVal = vec[i];
        }
    }
    
    void getMinMaxValues( const std::vector<VectorType> &vecs, RealType &minVal, RealType &maxVal ) const {
        minVal = 2.e+19; maxVal = -2.e+19;
        for( int i=0; i<vecs.size(); ++i ){
            for( int j=0; j<vecs[i].size(); ++j ){
                if( vecs[i][j] > maxVal ) maxVal = vecs[i][j];
                if( vecs[i][j] < minVal ) minVal = vecs[i][j];
            }
        }
    }
    
    void generateDataFile ( const VectorType &vec, const string &outputFileName ) const {
        std::ofstream outDataFile ( aol::strprintf ( "%s/%s.dat", _saveDirectory.c_str(), outputFileName.c_str() ).c_str()  );
        outDataFile << "x y" << endl;
        for( int i=0; i<vec.size(); ++i ) outDataFile << std::fixed << std::setprecision( 12 ) << i << " " << vec[i] << endl;
        outDataFile.close();

    }
    
    void generateDataFile ( const VectorType &t, const VectorType &vecs, const string &outputFileName ) const {
        std::ofstream outDataFile ( aol::strprintf ( "%s/%s.dat", _saveDirectory.c_str(), outputFileName.c_str() ).c_str()  );
        outDataFile << "x y" << endl;
        for( int j=0; j<t.size(); ++j ) outDataFile << std::fixed << std::setprecision( 12 ) << t[j] << " " << vecs[j] << endl;
        outDataFile.close();
    }
    
    void generateDataFiles ( const VectorType &t, const std::vector<VectorType> &vecs, const string &outputFileName ) const {
        for( int i=0; i<vecs.size(); ++i ){
            std::ofstream outDataFile ( aol::strprintf ( "%s/%s_%d.dat", _saveDirectory.c_str(), outputFileName.c_str(), i ).c_str()  );
            outDataFile << "x y" << endl;
            for( int j=0; j<t.size(); ++j ) outDataFile << std::fixed << std::setprecision( 12 ) << t[j] << " " << vecs[i][j] << endl;
            outDataFile.close();
        }
    }
    
    
    void plotCurve ( const VectorType &vec, const string &outputFileName, const string &title ) const {
        
        this->generateDataFile( vec, outputFileName );
        RealType minVal, maxVal; this->getMinMaxValues( vec, minVal, maxVal );
        
        string texFileName = aol::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outputFileName.c_str()  ).c_str ();
        std::ofstream out ( texFileName.c_str()  );
        
        _tikzHelper.generateIncludes( out, 10, 10, 0.5, 0.5, 0.5, 0.5 );
        _tikzHelper.generateBeginDocument( out );
        _tikzHelper.generateBeginTikzPicture( out, 0.9 );
        
        //plot
        out << "\\begin{axis}" << endl
            << " [ title=" << title.c_str() << "," << endl
            << "   width=0.8\\paperwidth," << endl
            << "   at={(0.0,0)}," << endl
            << "   xmin=0, xmax=" << vec.size() << "," << endl
            << "   ymin=0, ymax=" << maxVal << "," << endl
            << " ]" << endl;
        out << " \\addplot[thick] table {" << aol::strprintf ( "%s/%s.dat", _saveDirectory.c_str(), outputFileName.c_str() ).c_str() << "};" << endl;
        out << "\\end{axis}" << endl;
        // end plot
        
        _tikzHelper.generateEndTikzPicture( out );
        _tikzHelper.generateEndDocument( out );
        
        out.close();
        
        runPdflatex( outputFileName, _saveDirectory, true );
    }
    
    
    
   void plotCurve ( const VectorType &t, const VectorType &vec, const string &outputFileName, const string &title ) const {
        
        this->generateDataFile( t, vec, outputFileName );
        
        RealType min_x, max_x; this->getMinMaxValues( t, min_x, max_x );
        RealType minVal, maxVal; this->getMinMaxValues( vec, minVal, maxVal );
        
        string texFileName = aol::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outputFileName.c_str()  ).c_str ();
        std::ofstream out ( texFileName.c_str()  );
        
        _tikzHelper.generateIncludes( out, 10, 10, 0.5, 0.5, 0.5, 0.5 );
        _tikzHelper.generateBeginDocument( out );
        _tikzHelper.generateBeginTikzPicture( out, 0.9 );
        
        //plot
        out << "\\begin{axis}" << endl
            << " [ title=" << title.c_str() << "," << endl
            << "   width=0.8\\paperwidth," << endl
            << "   at={(0.0,0)}," << endl
            << "   xmin=" << min_x  << ", xmax=" << max_x << "," << endl
            << "   ymin=" << minVal << ", ymax=" << maxVal << "," << endl
            << "   legend style={at={(0.1,0.8)},anchor=south west}," << endl
            << " ]" << endl;
        out << " \\addplot[color = black] table {" << aol::strprintf ( "%s/%s.dat", _saveDirectory.c_str(), outputFileName.c_str() ).c_str() << "};" << endl;
        out << "\\end{axis}" << endl;
        // end plot
        
        _tikzHelper.generateEndTikzPicture( out );
        _tikzHelper.generateEndDocument( out );
        
        out.close();
        
        runPdflatex( outputFileName, _saveDirectory, true );
    }
    
    
    void plotMutlipleCurves ( const VectorType &t, const std::vector<VectorType> &vecs, const string &outputFileName, const string &title,
        std::vector<string> &colors, std::vector<string> &describtions) const {
        
        this->generateDataFiles( t, vecs, outputFileName );
        
        RealType min_x, max_x; this->getMinMaxValues( t, min_x, max_x );
        RealType minVal, maxVal; this->getMinMaxValues( vecs, minVal, maxVal );
        
        string texFileName = aol::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outputFileName.c_str()  ).c_str ();
        std::ofstream out ( texFileName.c_str()  );
        
        _tikzHelper.generateIncludes( out, 10, 10, 0.5, 0.5, 0.5, 0.5 );
        _tikzHelper.generateBeginDocument( out );
        _tikzHelper.generateBeginTikzPicture( out, 0.9 );
        
        //plot
        out << "\\begin{axis}" << endl
            << " [ title=" << title.c_str() << "," << endl
            << "   width=0.8\\paperwidth," << endl
            << "   at={(0.0,0)}," << endl
            << "   xmin=" << min_x  << ", xmax=" << max_x << "," << endl
            << "   ymin=" << minVal << ", ymax=" << maxVal << "," << endl
            << "   legend style={at={(0.1,0.8)},anchor=south west}," << endl
            << " ]" << endl;
        for( int i=0; i<vecs.size(); ++i ){
            out << " \\addplot[color = " << colors[i].c_str() << "] table {" << aol::strprintf ( "%s/%s_%d.dat", _saveDirectory.c_str(), outputFileName.c_str(), i ).c_str() << "};" << endl;
        }
        out << "  \\legend{";
        for( int i=0; i<vecs.size(); ++i ){
             out << describtions[i].c_str();
             if( i < vecs.size() - 1 ) out << ", ";
             else out << "}" << endl;
        }
        out << "\\end{axis}" << endl;
        // end plot
        
        _tikzHelper.generateEndTikzPicture( out );
        _tikzHelper.generateEndDocument( out );
        
        out.close();
        
        runPdflatex( outputFileName, _saveDirectory, true );
    }
    
    
    void plotMutlipleCurves ( const std::vector<VectorType> &t, const std::vector<VectorType> &vecs, const string &outputFileName, const string &title,
        std::vector<string> &colors, std::vector<string> &describtions) const {

        for( int i=0; i<vecs.size(); ++i ){
            std::ofstream outDataFile ( aol::strprintf ( "%s/%s_%d.dat", _saveDirectory.c_str(), outputFileName.c_str(), i ).c_str()  );
            outDataFile << "x y" << endl;
            for( int j=0; j<t[i].size(); ++j ) outDataFile << std::fixed << std::setprecision( 12 ) << t[i][j] << " " << vecs[i][j] << endl;
            outDataFile.close();
        }
        
        RealType min_x, max_x; this->getMinMaxValues( t, min_x, max_x );
        RealType minVal, maxVal; this->getMinMaxValues( vecs, minVal, maxVal );
        
        string texFileName = aol::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outputFileName.c_str()  ).c_str ();
        std::ofstream out ( texFileName.c_str()  );
        
        _tikzHelper.generateIncludes( out, 10, 10, 0.5, 0.5, 0.5, 0.5 );
        _tikzHelper.generateBeginDocument( out );
        _tikzHelper.generateBeginTikzPicture( out, 0.9 );
        
        //plot
        out << "\\begin{axis}" << endl
            << " [ title=" << title.c_str() << "," << endl
            << "   width=0.8\\paperwidth," << endl
            << "   at={(0.0,0)}," << endl
            << "   xmin=" << min_x  << ", xmax=" << max_x << "," << endl
            << "   ymin=" << minVal << ", ymax=" << maxVal << "," << endl
            << "   legend style={at={(0.1,0.8)},anchor=south west}," << endl
            << " ]" << endl;
        for( int i=0; i<vecs.size(); ++i ){
            out << " \\addplot[color = " << colors[i].c_str() << "] table {" << aol::strprintf ( "%s/%s_%d.dat", _saveDirectory.c_str(), outputFileName.c_str(), i ).c_str() << "};" << endl;
        }
        out << "  \\legend{";
        for( int i=0; i<vecs.size(); ++i ){
             out << describtions[i].c_str();
             if( i < vecs.size() - 1 ) out << ", ";
             else out << "}" << endl;
        }
        out << "\\end{axis}" << endl;
        // end plot
        
        _tikzHelper.generateEndTikzPicture( out );
        _tikzHelper.generateEndDocument( out );
        
        out.close();
        
        runPdflatex( outputFileName, _saveDirectory, true );
    }
    
    void plotLogLogCurve ( const VectorType &vec, const string &outputFileName, const string &title ) const {
        
        this->generateDataFile( vec, outputFileName );
        
        string texFileName = aol::strprintf ( "%s/%s.tex", _saveDirectory.c_str(), outputFileName.c_str()  ).c_str ();
        std::ofstream out ( texFileName.c_str()  );
        
        _tikzHelper.generateIncludes( out, 10, 10, 0.5, 0.5, 0.5, 0.5 );
        _tikzHelper.generateBeginDocument( out );
        _tikzHelper.generateBeginTikzPicture( out, 0.9 );  
        
        //loglogPlot
        out << "\\begin{loglogaxis}" << endl
            << " [ xlabel = \\textsc{iteration}," << endl
            << "   %ylabel = $$" << endl
            << " ]" << endl;
        out << " \\addplot[thick] table {" << aol::strprintf ( "%s/%s.dat", _saveDirectory.c_str(), outputFileName.c_str() ).c_str() << "};" << endl;
        out << " %\\legend{$$}" << endl;
        out << "\\end{loglogaxis}" << endl;
        //end loglogPlot
        
        _tikzHelper.generateEndTikzPicture( out );
        _tikzHelper.generateEndDocument( out );
        
        out.close();
        
        runPdflatex( outputFileName, _saveDirectory, true );
    }

};


#endif //__TIKZPLOTTER_H