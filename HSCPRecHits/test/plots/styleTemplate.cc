#include "TStyle.h"
TStyle *styleTemplate;

// tdrGrid: Turns the grid lines on (true) or off (false)

void tdrGrid(bool gridOn) {
  styleTemplate->SetPadGridX(gridOn);
  styleTemplate->SetPadGridY(gridOn);
}

// fixOverlay: Redraws the axis

void fixOverlay() {
  gPad->RedrawAxis();
}

void setStyleTemplate() {
  styleTemplate = new TStyle("styleTemplate","Style for P-TDR");

// For the canvas:
  styleTemplate->SetCanvasBorderMode(0);
  styleTemplate->SetCanvasColor(kWhite);
  styleTemplate->SetCanvasDefH(600); //Height of canvas
  styleTemplate->SetCanvasDefW(600); //Width of canvas
  styleTemplate->SetCanvasDefX(0);   //POsition on screen
  styleTemplate->SetCanvasDefY(0);

// For the Pad:
  styleTemplate->SetPadBorderMode(0);
  // styleTemplate->SetPadBorderSize(Width_t size = 1);
  styleTemplate->SetPadColor(kWhite);
  styleTemplate->SetPadGridX(false);
  styleTemplate->SetPadGridY(false);
  styleTemplate->SetGridColor(0);
  styleTemplate->SetGridStyle(3);
  styleTemplate->SetGridWidth(1);

// For the frame:
  styleTemplate->SetFrameBorderMode(0);
  styleTemplate->SetFrameBorderSize(1);
  styleTemplate->SetFrameFillColor(0);
  styleTemplate->SetFrameFillStyle(0);
  styleTemplate->SetFrameLineColor(1);
  styleTemplate->SetFrameLineStyle(1);
  styleTemplate->SetFrameLineWidth(1);
  
// For the histo:
  // styleTemplate->SetHistFillColor(1);
  // styleTemplate->SetHistFillStyle(0);
  styleTemplate->SetHistLineColor(1);
  styleTemplate->SetHistLineStyle(0);
  styleTemplate->SetHistLineWidth(1);
  // styleTemplate->SetLegoInnerR(Float_t rad = 0.5);
  // styleTemplate->SetNumberContours(Int_t number = 20);

  styleTemplate->SetEndErrorSize(2);
  // styleTemplate->SetErrorMarker(20);
  //styleTemplate->SetErrorX(0.);
  
  styleTemplate->SetMarkerStyle(20);
  
//For the fit/function:
  styleTemplate->SetOptFit(1);
  styleTemplate->SetFitFormat("5.4g");
  styleTemplate->SetFuncColor(2);
  styleTemplate->SetFuncStyle(1);
  styleTemplate->SetFuncWidth(1);

//For the date:
  styleTemplate->SetOptDate(0);
  // styleTemplate->SetDateX(Float_t x = 0.01);
  // styleTemplate->SetDateY(Float_t y = 0.01);

// For the statistics box:
  styleTemplate->SetOptFile(0);
  styleTemplate->SetOptStat(0); // To display Entries, Mean, RMS, Underflow and OverFlow styleTemplate->SetOptStat("nemruo");
  styleTemplate->SetStatColor(kWhite);
  styleTemplate->SetStatFont(42);
  styleTemplate->SetStatFontSize(0.025);
  styleTemplate->SetStatTextColor(1);
  styleTemplate->SetStatFormat("6.4g");
  styleTemplate->SetStatBorderSize(1);
  styleTemplate->SetStatH(0.1);
  styleTemplate->SetStatW(0.15);
  // styleTemplate->SetStatStyle(Style_t style = 1001);
  // styleTemplate->SetStatX(Float_t x = 0);
  // styleTemplate->SetStatY(Float_t y = 0);

// Margins:
//  styleTemplate->SetPadTopMargin(0.05);
//  styleTemplate->SetPadBottomMargin(0.13);
//  styleTemplate->SetPadLeftMargin(0.13);
//  styleTemplate->SetPadRightMargin(0.05);

// For the Global title:

  styleTemplate->SetOptTitle(0);
  styleTemplate->SetTitleFont(42);
  styleTemplate->SetTitleColor(1);
  styleTemplate->SetTitleTextColor(1);
  styleTemplate->SetTitleFillColor(10);
  styleTemplate->SetTitleFontSize(0.05);
  // styleTemplate->SetTitleH(0); // Set the height of the title box
  // styleTemplate->SetTitleW(0); // Set the width of the title box
  // styleTemplate->SetTitleX(0); // Set the position of the title box
  // styleTemplate->SetTitleY(0.985); // Set the position of the title box
  // styleTemplate->SetTitleStyle(Style_t style = 1001);
  // styleTemplate->SetTitleBorderSize(2);

// For the axis titles:

  styleTemplate->SetTitleColor(1, "XYZ");
  styleTemplate->SetTitleFont(42, "XYZ");
  styleTemplate->SetTitleSize(0.06, "XYZ");
  // styleTemplate->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // styleTemplate->SetTitleYSize(Float_t size = 0.02);
  styleTemplate->SetTitleXOffset(1.0);
  styleTemplate->SetTitleYOffset(1.25);
  styleTemplate->SetTitleXSize(0.04);
  styleTemplate->SetTitleYSize(0.04);
  // styleTemplate->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  styleTemplate->SetLabelColor(1, "XYZ");
  styleTemplate->SetLabelFont(42, "XYZ");
  styleTemplate->SetLabelOffset(0.007, "XYZ");
  styleTemplate->SetLabelSize(0.04, "XYZ");

// For the axis:
  //TGaxis::SetMaxDigits(3);
  styleTemplate->SetAxisColor(1, "XYZ");
  //styleTemplate->SetStripDecimals(kTRUE);
  styleTemplate->SetTickLength(0.03, "XYZ");
  styleTemplate->SetNdivisions(510, "XYZ");
  styleTemplate->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  styleTemplate->SetPadTickY(1);

// Change for log plots:
  styleTemplate->SetOptLogx(0);
  styleTemplate->SetOptLogy(0);
  styleTemplate->SetOptLogz(0);

// Postscript options:
  styleTemplate->SetPaperSize(20.,20.);
  // styleTemplate->SetLineScalePS(Float_t scale = 3);
  // styleTemplate->SetLineStyleString(Int_t i, const char* text);
  // styleTemplate->SetHeaderPS(const char* header);
  // styleTemplate->SetTitlePS(const char* pstitle);

  // styleTemplate->SetBarOffset(Float_t baroff = 0.5);
  // styleTemplate->SetBarWidth(Float_t barwidth = 0.5);
  // styleTemplate->SetPaintTextFormat(const char* format = "g");
  // styleTemplate->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // styleTemplate->SetTimeOffset(Double_t toffset);
  // styleTemplate->SetHistMinimumZero(kTRUE);

  styleTemplate->SetHatchesLineWidth(5);
  styleTemplate->SetHatchesSpacing(0.05);

  
  styleTemplate->cd();

}
