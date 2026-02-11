CustomComponents <- tags$script(HTML(sprintf("(function() {
  const React = jsmodule['react'];
  const ReactDOMServer = jsmodule['react-dom'];
  const Shiny = jsmodule['@/shiny'];
  const CustomComponents = jsmodule['CustomComponents'] ??= {};

  CustomComponents.TaxoniumComponent = function(propz) {
  
  //console.log('META TEST')
  //console.log(propz)
  
  propz.treeData[0].metadata = propz.meta[0]; //enable metada annotation
  
  /*
  Colour palette for H x
  1.	[1, 115, 178],
	2.	[222, 143, 5],
	3.	[2, 158, 115],
	4.	[213, 94, 0],
	5.	[204, 120, 188],
	6.	[202, 145, 97],
	7.	[251, 175, 228],
	8.	[148, 148, 148],
	9.	[236, 225, 51],
	
	10.	[86, 180, 233],
	11.	[0, 0, 255], #(Blue)
	12.	[255, 0, 0], #(Red)
	13.	[0, 255, 0], #(Green)
	14.	[255, 255, 0], #(Yellow)
	
	15.	[255, 0, 255], #(Magenta)
	16.	[0, 255, 255], #(Cyan)
	17.	[128, 0, 0], #(Dark Red)
	18.	[0, 128, 0], #(Dark Green)
	
	19.	[0, 0, 128], #(Dark Blue)
	20.	[128, 128, 0], #(Brown)
	21.	[0, 128, 128], #(Teal)
	22.	[128, 0, 128], #(Purple)
	23.	[255, 128, 0], #(Orange)
	24.	[0, 128, 255] #(Azure)
	
	Colour palette for N subtype
Here are the RGB values for the colorblind-safe palette for 13 categories, printed in the original format:

1. **Category 1**: RGB(1, 115, 178)
2. **Category 2**: RGB(222, 143, 5)
3. **Category 3**: RGB(2, 158, 115)
4. **Category 4**: RGB(213, 94, 0)
5. **Category 5**: RGB(204, 120, 188)

6. **Category 6**: RGB(202, 145, 97)
7. **Category 7**: RGB(251, 175, 228)
8. **Category 8**: RGB(148, 148, 148)
9. **Category 9**: RGB(0, 114, 178)
10. **Category 10**: RGB(230, 159, 0)

11. **Category 11**: RGB(86, 180, 233)
12. **Category 12**: RGB(240, 228, 66)
13. **Category 13**: RGB(102, 166, 30)

Colour palette for Host group - Brewer palette Dark2 as RGB values
1.	#1B9E77: RGB(27, 158, 119)
2.	#D95F02: RGB(217, 95, 2)
3.	#7570B3: RGB(117, 112, 179)
4.	#E7298A: RGB(231, 41, 138)
5.	#66A61E: RGB(102, 166, 30)

Colour palette for host order
  1.	Category 1: [1, 115, 178],
	2.	Category 2: [222, 143, 5],
	3.	Category 3: [2, 158, 115],
	4.	Category 4: [213, 94, 0],
	5.	Category 5: [204, 120, 188],
	6.	Category 6: [202, 145, 97],
	7.	Category 7: [251, 175, 228],
	8.	Category 8: [148, 148, 148],
	9.	Category 9: [1, 138, 213],
	10.	Category 10: [0, 92, 142],
	11.	Category 11: [1, 103, 178],
	12.	Category 12: [0, 126, 178],
	13.	Category 13: [255, 171, 6],
	14.	Category 14: [177, 114, 4],
  15.	Category 15: [244, 128, 5],
	16.	Category 16: [199, 157, 5],
	17.	Category 17: [2, 189, 138],
	18.	Category 18: [1, 126, 92],
	19.	Category 19: [2, 142, 115],
	20.	Category 20: [1, 173, 115],
	21.	Category 21: [255, 112, 0],
	22.	Category 22: [170, 75, 0],
	23.	Category 23: [234, 84, 0],
	24.	Category 24: [191, 103, 0],
	25.	Category 25: [244, 144, 225],
	26.	Category 26: [163, 96, 150],
	27.	Category 27: [224, 108, 188],
	28.	Category 28: [183, 132, 188],
	
	*/
  // Amino acid colour palette - standardized chemically-indexed colors
  const aaPalette = JSON.parse('%s');
  
  const config = {'colorMapping':{
  // H subtype colour palette
  'H1':[1, 115, 178],
  'H3':[222, 143, 5],
  'H9':[2, 158, 115],
  'H5':[213, 94, 0],
  'H6':[204, 120, 188],
  
  'H7':[202, 145, 97],
  'NA':[251, 175, 228],
  'H4':[148, 148, 148],
  'H10':[236, 225, 51],
  
  'H11':[86, 180, 233],
  'H13':[0, 0, 255],
  'H2':[255, 0, 0],
  'H16':[0, 255, 0],
  'H12':[255, 255, 0],
  
  'H8':[255, 0, 255],
  'H18':[0, 255, 255],
  'H14':[128, 0, 0],
  'H15':[0, 128, 0],
  
  'H17':[0, 0, 128],
  '':[255, 255, 255],
  'Hx':[0, 128, 128],
  'H19':[128, 0, 128],
  'H1n2':[255, 128, 0],
  'unknown':[0, 128, 255],
  
  // N subtype colour palette
  'N2':[1, 115, 178],
  'N1':[222, 143, 5],
  'N6':[2, 158, 115],
  'N8':[213, 94, 0],
  'NA':[204, 120, 188],
  
  'N3':[202, 145, 97],
  'N9':[251, 175, 228],
  'N7':[148, 148, 148],
  'N5':[0, 114, 178],
  'N4':[230, 159, 0],
  
  'N11':[86, 180, 233],
  'N10':[240, 228, 66],
  'Nx':[102, 166, 30],
  
  // Host group - Brewer colour palette Set2
  'Birds':[102, 194, 165],
  'Human':[252, 141, 98],
  'Other Mammals':[141, 160, 203],
  'Environment':[231, 138, 195],
  'NA':[166, 216, 84],
  'Unknown':[255, 217, 47],
  //'Unknown':[229, 196, 148],
  //'Unknown':[179, 179, 179],
  
  // Host order - extended Seaborn colour blind palette
  'Artiodactyla':[1, 115, 178],
  'Anseriformes':[222, 143, 5],  
  'Galliformes':[2, 158, 115],
  //'Unknown':[213, 94, 0],
  
  'Unknown avian':[204, 120, 188], 
  'Charadriiformes':[202, 145, 97],
  'Primates':[251, 175, 228],
  'Environment':[148, 148, 148],
  
  'Chiroptera':[1, 138, 213],
  'Carnivora':[0, 92, 142],
  'Perissodactyla':[1, 103, 178],
  'Columbiformes':[0, 126, 178],
  
  'Struthioniformes':[255, 171, 6],
  'Gruiformes':[177, 114, 4],
  'Passeriformes':[244, 128, 5],
  'Pelecaniformes':[199, 157, 5],
  
  'Sphenisciformes':[2, 189, 138],
  'Accipitriformes':[1, 126, 92],
  'Casuariiformes':[2, 142, 115],
  'Phoenicopteriformes':[1, 173, 115],
  
  'Lagomorpha':[255, 112, 0],
  'Podicipediformes':[170, 75, 0],
  'Suliformes':[234, 84, 0],
  'Otidiformes':[191, 103, 0],
  
  'Procellariiformes':[244, 144, 225],
  'Psittaciformes':[163, 96, 150],
  'Strigiformes':[224, 108, 188], 
  'Tinamiformes':[183, 132, 188],
  
  ...aaPalette
  }};
  
  // pass key to rerender child component when tree file or selected rows updated
   return React.createElement(Taxonium, { 
   sourceData: propz.treeData[0], 
   configDict: config,
   key: propz.treeData[0].filename.concat(propz.treeData[0].metadata.rows) 
   })
  };
  
})()", get_aa_palette_js()))) # end script

TaxoniumComponent <- function(...) {
  shiny.react::reactElement(
    module = "CustomComponents",
    name = "TaxoniumComponent",
    props = shiny.react::asProps(...),
  )
}
