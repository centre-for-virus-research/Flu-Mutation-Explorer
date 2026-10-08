# sprintf caps a format string at 8192 characters and this script is well past
# that, so the palette is pasted in rather than formatted in
CustomComponents <- tags$script(HTML(paste0("(function() {
  const React = jsmodule['react'];
  const ReactDOMServer = jsmodule['react-dom'];
  const Shiny = jsmodule['@/shiny'];
  const CustomComponents = jsmodule['CustomComponents'] ??= {};

  CustomComponents.TaxoniumComponent = function(propz) {
  
  propz.treeData[0].metadata = propz.meta[0]; // enable metadata annotation

  // The colour mappings below are colourblind-safe palettes: H and N subtype and host
  // order from the extended Seaborn colorblind palette, host group from Brewer Set2.

  // Amino acid colour palette - standardized chemically-indexed colors
  const aaPalette = JSON.parse('", get_aa_palette_js(), "');
  
  // Fields shown on the tips but not offered as a colour scheme. The cluster
  // composition strings are free text with a distinct value per cluster, so
  // colouring by them produces one colour per tip and no usable legend.
  const notColourable = ['cluster_members', 'cluster_hosts'];
  
  // keys_to_display (the pop-up) and colorBy.colorByOptions (the menu) are
  // separate config fields, and useConfig applies configDict last, so overriding
  // colorBy here drops fields from the menu while leaving the pop-up untouched.
  //
  // The columns differ between the default tree, a position search and a BLAST
  // hit selection, so the list is read from the metadata actually being rendered.
  // Header parsing mirrors Taxonium's own CSV handling: split on commas and strip
  // quotes (values are comma-free upstream). Column one is the tip name.
  const metaHeader = (propz.meta[0].data || '').split('\\n')[0].split(',')
    .map(function(h) { return h.replace(/\"/g, '').trim(); });
  
  const colourByOptions = metaHeader.slice(1)
    .filter(function(h) { return h && notColourable.indexOf(h) === -1; })
    .map(function(h) { return 'meta_' + h; });
  colourByOptions.push('None');
  // Strains of the BLAST hits the user selected, fed into Taxonium's own name
  // search so the tips are found and zoomed to instead of being left for the
  // user to hunt down. Shiny sends a character vector, which arrives as an
  // array, but a one-element vector can arrive as a bare string.
  const rawSearch = propz.search;
  const hitNames =
    (rawSearch == null ? [] : (Array.isArray(rawSearch) ? rawSearch : [rawSearch]))
      .filter(function(name) { return typeof name === 'string' && name.trim() !== ''; });
  const hasHits = hitNames.length > 0;
  
  // Taxonium keeps its search state in query.srch, a JSON string holding an
  // array of search specs. text_per_line matches node.name exactly - lowercased
  // and trimmed, one name per line - which is what a multi-row hits selection
  // needs, and the tree tip labels are the strain names the hits table shows.
  //
  // The hit search takes slot 0 - the slot the user would otherwise paste into.
  // Result circles are coloured by position in this array against a fixed palette
  // [[255,0,0],[0,0,255],[0,255,0],...], so slot 0 is the red one and anywhere
  // else comes out blue or green.
  const blankSpec = {key: 'taxonium_default', type: 'name', method: 'text_match',
                     text: '', gene: 'S', position: 484, new_residue: 'any', min_tips: 0};
  const hitsSpec  = {key: 'blast_hits', type: 'name', method: 'text_per_line',
                     text: hitNames.join('\\n'), gene: 'S', position: 484,
                     new_residue: 'any', min_tips: 0};
  
  const initialQuery = {
    srch: JSON.stringify(hasHits ? [hitsSpec] : [blankSpec]),
    enabled: JSON.stringify(hasHits ? {blast_hits: true} : {taxonium_default: true}),
    backend: '',
    xType: 'x_dist',
    mutationTypesEnabled: JSON.stringify({aa: true, nt: false}),
    treenomeEnabled: false
  };
  
  // Zoom to the hit search. The index goes in as the string '0', not the number:
  // this build gates the zoom on `query.zoomToSearch ? {index: ...} : null`, so a
  // numeric 0 is falsy and would never zoom, while '0' is truthy and still finds
  // searchSpec[0] because the array lookup coerces it back.
  //
  // The zoom is what keeps a place in the tree across a remount: useConfig runs
  // once per mount and calls onViewStateChange with the tree's own initial_y and
  // initial_zoom, so every mount otherwise resets the view to the whole tree. The
  // remount cannot be avoided either - the node re-query depends on config, which
  // only useConfig sets, so metadata arriving without a mount is uploaded to the
  // worker and never queried, leaving the tree on its old colours.
  //
  // Zooming narrows the colour key as a side effect, since Taxonium builds the key
  // from the nodes in the viewport. Zooming back out restores it; landing on the
  // hit was judged worth more than the wider key.
  if(hasHits) { initialQuery.zoomToSearch = '0'; }
  
  // query and updateQuery have to be supplied together: this build falls back to
  // its own state only when both are absent, so passing query alone would leave
  // the search panel unable to edit itself. Taxonium also calls updateQuery to
  // clear zoomToSearch once it has zoomed, so this has to be real state.
  const [query, setQuery] = React.useState(initialQuery);
  const updateQuery = React.useCallback(function(patch) {
    setQuery(function(previous) { return Object.assign({}, previous, patch); });
  }, []);
  
  // The key below remounts the Taxonium child, but Shiny re-renders this wrapper in
  // place, and a useState initialiser only runs on the first mount. So the query is
  // re-seeded by hand whenever the child is about to remount, during render rather
  // than from an effect, which puts the new query and the new key in the same render
  // - Taxonium reads zoomToSearch only while mounting.
  //
  // The seed is everything the key is built from, so the two move together. Two
  // pieces of state ride on it:
  //
  //  - zoomToSearch, which Taxonium clears once it has zoomed. Left cleared, a
  //    remount falls back to the whole tree and a position search throws away the
  //    zoom to the hit.
  //  - color, the colour-by field, which outranks defaultColorByField in useColorBy
  //    and would otherwise survive every remount, so a position search would keep
  //    colouring by subtype rather than by the amino acids just loaded. Dropped when
  //    the columns change, since the field it names may no longer exist, and kept
  //    otherwise so a colour picked by hand survives selecting another hit.
  const remountSeed = propz.treeData[0].filename + '|' + propz.meta[0].rows +
                      '|' + hitNames.join('|'); // the key, verbatim
  const colourSeed  = colourByOptions.join('|');
  
  const [seen, setSeen] = React.useState({remount: remountSeed, colours: colourSeed});
  
  if(seen.remount !== remountSeed || seen.colours !== colourSeed) {
    const sameColumns = seen.colours === colourSeed;
    setSeen({remount: remountSeed, colours: colourSeed});
    setQuery(function(previous) {
      const next = Object.assign({}, initialQuery);
      if(sameColumns && previous.color !== undefined) { next.color = previous.color; }
      return next;
    });
  }
  
  
  const config = {
  'colorBy': {'colorByOptions': colourByOptions},
  // Taxonium would default to the first generated option, which may be one that
  // was just removed above
  'defaultColorByField': colourByOptions[0],
  'colorMapping':{
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
  
   return React.createElement(Taxonium, {
   sourceData: propz.treeData[0], 
   configDict: config,
   query: query,
   updateQuery: updateQuery,
   // The key remounts the child when the tree, the metadata or the hits change.
   // The hit names are part of it because a second BLAST run can select the same
   // row numbers as the first, leaving metadata.rows unchanged while the strains
   // behind them differ - without a remount the new search would never apply.
   key: propz.treeData[0].filename
          .concat(propz.treeData[0].metadata.rows)
          .concat(hitNames.join('|'))
   })
  };
  
})()"))) # end script

TaxoniumComponent <- function(...) {
  shiny.react::reactElement(
    module = "CustomComponents",
    name = "TaxoniumComponent",
    props = shiny.react::asProps(...),
  )
}
