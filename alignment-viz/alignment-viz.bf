RequireVersion ("2.5.40");

LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/convenience/regexp.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");
LoadFunctionLibrary     ("lib/post-script.bf");


utility.SetEnvVariable  ("NORMALIZE_SEQUENCE_NAMES", TRUE);
utility.SetEnvVariable  ("ACCEPT_ROOTED_TREES", TRUE);


al_viz.analysis_description = {terms.io.info :
                            "
                                Load a sequence alignment (nucleotide or A/A), and generate a PostScript visualization of it.
                            ",
                            terms.io.version :          "0.01",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "A sequence alignment"
                          };

io.DisplayAnalysisBanner (al_viz.analysis_description);

KeywordArgument ("alignment", "The alignment to visualize");
av_viz.alignment = alignments.ReadNucleotideDataSet ("av_viz.raw_data", None);

io.ReportProgressMessage ("Data QC", "Loaded `av_viz.alignment[terms.data.sequences]` sequences on `av_viz.alignment[terms.data.sites]` sites from **`av_viz.alignment[terms.data.file]`**");

KeywordArgument ("site-filter", "Use this string to filter sites", "None");
av_viz.site_filter = io.PromptUserForString ("Use this string to filter sites");
if (av_viz.site_filter == "None") {
    av_viz.site_filter = "";
}

KeywordArgument ("sequence-filter", "Use this string to filter sequences", "None");
av_viz.seq_filter = io.PromptUserForString ("Use this string to filter sequences");

if (av_viz.seq_filter == "None") {
    av_viz.seq_filter = "";
}

KeywordArgument ("sequence-regexp", "Use this regexp to filter sequences by name", "None");
av_viz.seq_regexp = io.PromptUserForString ("Use this regexp to filter sequences by name");



if (av_viz.seq_regexp != "None") {
    function av_viz.filter_regexp (id, data) {
        return None != regexp.Find (id, av_viz.seq_regexp);
    }

    DataSetFilter av_viz.filter = CreateFilter (av_viz.raw_data,1, av_viz.site_filter, "av_viz.filter_regexp");
} else {
    DataSetFilter av_viz.filter = CreateFilter (av_viz.raw_data,1, av_viz.site_filter, av_viz.seq_filter);
}

io.ReportProgressMessage ("Data QC", "Filtering retained `av_viz.filter.sites` sites and `av_viz.filter.species` sequences");



KeywordArgument ("output", "Write the PostScript to", av_viz.alignment[terms.data.file] + ".ps");
av_viz.output = io.PromptUserForFilePath ("Write the PostScript to");

KeywordArgument ("font-size", "Baseline font size [4-36]", 14);
av_viz.font_size = io.PromptUser ("Baseline font size", 14,4,36,1);

KeywordArgument ("ticks", "Display site indices every X sites [3-20]", 5);
av_viz.ticks = io.PromptUser ("Display site indices every X sites", 5,3,20,1);



_charColorList = {};
_charColorList["A"] = {{31, 119, 180}}*(1/255);
_charColorList["C"] = {{255, 127, 14}}*(1/255);
_charColorList["D"] = {{102,0,153}}*(1/255);
_charColorList["E"] = {{153,0,153}}*(1/255);
_charColorList["F"] = {{255,51,51}}*(1/255);
_charColorList["G"] = {{44, 160, 44}}*(1/255);
_charColorList["H"] = {{255,204,0}}*(1/255);
_charColorList["I"] = {{153,255,153}}*(1/255);
_charColorList["K"] = {{102,153,51}}*(1/255);
_charColorList["L"] = {{102,204,102}}*(1/255);
_charColorList["M"] = {{0,51,153}}*(1/255);
_charColorList["N"] = {{255,102,0}}*(1/255);
_charColorList["P"] = {{153,153,51}}*(1/255);
_charColorList["Q"] = {{255,153,0}}*(1/255);
_charColorList["R"] = {{153,51,102}}*(1/255);
_charColorList["S"] = {{153,153,153}}*(1/255);
_charColorList["T"] = {{214, 39, 40}}*(1/255);
_charColorList["V"] = {{153,51,51}}*(1/255);
_charColorList["W"] = {{204,51,51}}*(1/255);
_charColorList["Y"] = {{102,153,102}}*(1/255);

for (letter, color; in; _charColorList) {
    _charColorList[letter] = Join (" ", color);
}

_font_size  = av_viz.font_size;
_char_space_mult = 1 + 0.1*((Log (av_viz.filter.sites)/Log(2) + 1.5)$1);
_char_space = (_font_size*_char_space_mult);
_page_w	    = 612;
_page_h     = 792;
_atom		= av_viz.ticks;

fprintf (av_viz.output,CLEAR_FILE,_HYPSPageHeader (_page_w,_page_h,"Character Plot for "+ av_viz.alignment[terms.data.file]),
							_HYPSSetFont ("Courier",_font_size),
							"/dmx 6 array currentmatrix def\n",
							"/sshift fs 2 idiv fs add def\n",
							"/setb {0 0 0 setrgbcolor} def\n",
							"/flbx {setrgbcolor newpath 2 copy moveto 4 -1 roll exch 2 copy lineto exch 4 -1 roll 2 copy lineto 4 -1 roll exch lineto pop pop closepath fill setb} def\n",
							"/drawletter {translate scale newpath 0 0 moveto false charpath fill dmx setmatrix translate 0.4 0.4 scale newpath sshift 0 moveto  false charpath fill dmx setmatrix} def\n"
);



GetDataInfo (_charHandles, av_viz.filter, "CHARACTERS");
_ccount  = Columns (_charHandles);
_unit    = {1,_ccount}["1"];
GetDataInfo (_dupInfo, av_viz.filter);
_result_cache = {};
_maxD		  = 0;

_char_per_line  = _page_w / _char_space $ _atom * _atom;
if (_char_per_line == 0) {
	io.ReportAnExecutionError ("\nERROR: At least " + _atom + " characters must fit in a line; reduce font size");
}
else {
	fprintf (stdout, "\n", _char_per_line, " characters per line\n");
}

_dbyLine = {};

for (_idx = 0; _idx < av_viz.filter.sites; _idx += 1) {
	_siteInfo = {_ccount, 2};
	_cCounter = {_ccount, 1};
	for (_sidx = 0; _sidx < av_viz.filter.species; _sidx += 1) {
		GetDataInfo (_thisChar, av_viz.filter, _sidx, _dupInfo[_idx]);
		/* don't count gaps */
		if (Abs (_thisChar) != Sqrt (_ccount)) {
			_cCounter = _cCounter + _thisChar*(1/(_unit*_thisChar)[0]);
		}
	}
	_siteInfo = _siteInfo ["_MATRIX_ELEMENT_ROW_ * _MATRIX_ELEMENT_COLUMN_ + (1-_MATRIX_ELEMENT_COLUMN_)*_cCounter[_MATRIX_ELEMENT_ROW_]"]%0;
	for (_sidx = _ccount-1; _sidx >= 0; _sidx = _sidx - 1) {
		if (_siteInfo[_sidx][0] == 0) {
			break;
		}
	}
	
	_result_cache[_idx] = _siteInfo[{{_sidx+1,0}}][{{_ccount-1,1}}];
	_sidx = Rows (_result_cache[_idx]);
	
	if (_sidx > _maxD) {
		_maxD = _sidx;
	}

	if ((_idx + 1)%_char_per_line==0) {
		_dbyLine [Abs(_dbyLine)] = _maxD;
		_maxD = 0;
	}	
}

_current_x 	  = 0;
_current_y	  = _page_h-2*_font_size;
_cl			  = 0;

for (_idx = 0; _idx < av_viz.filter.sites; _idx += 1) {
	_cCounter = _result_cache[_idx];
	
	if (Rows(_cCounter)) {
		for (_sidx = Rows(_cCounter)-1; _sidx >= 0; _sidx = _sidx - 1) {
			_siteInfo = _current_y-_font_size*(Rows(_cCounter)-1-_sidx);
			_my_c     = _charHandles[_cCounter[_sidx][1]];
			if (Abs(_charColorList[_my_c])){
				fprintf (av_viz.output, _charColorList[_my_c], " setrgbcolor\n");
			}
			
			_char_count = _cCounter[_sidx][0];
			if (_char_count $ 1 != _char_count) {
			    _char_count = Format (_char_count, -1, 2);
			}
			
			fprintf (av_viz.output, "(", _char_count ,") ", _current_x, " ",_siteInfo," (",_my_c,") 1 1 ", _current_x, " ", _siteInfo, " drawletter\n");
			if (Abs(_charColorList[_my_c])) {
				fprintf (av_viz.output, "setb\n");
			}
		}
	}
	else {
		fprintf (av_viz.output, "( ) ", _current_x, " ",_current_y," (-) 1 1 ", _current_x, " ", _current_y, " drawletter\n");
	}
	_current_x = _current_x + _char_space;

	if ((_idx + 1)%_char_per_line==0 || _idx == av_viz.filter.sites-1) {
		_current_y = _current_y + _font_size;
		if (_idx == av_viz.filter.sites-1)
		{
			if (av_viz.filter.sites % _char_per_line == 0) {
				_startx = _char_per_line;
			}
			else {
				_startx = av_viz.filter.sites%_char_per_line$_atom*_atom;
			}
			_idx2 = av_viz.filter.sites-av_viz.filter.sites%_char_per_line;
		}
		else {
			_idx2 = _idx - _char_per_line + 1;
			_startx = _char_per_line;
		}
		fprintf (av_viz.output, "0 ", _current_y + _font_size * 4$5, " ", (_char_space+1) * _char_per_line , " ", _current_y - _font_size$4, " 0.9 0.9 0.9 flbx\n");
		for (_idx3 = _startx; _idx3 > 0; _idx3 = _idx3 - _atom) {
			fprintf (av_viz.output, "( ) 0 0 (",_idx2+_idx3,") 0.9 0.9 ", (_idx3-1) * _char_space, " ", _current_y, " drawletter\n");
		}
		_current_x = 0; 
		_current_y = _current_y - (2+_dbyLine[_cl])*_font_size;
		_cl = _cl+1;
		if (_current_y - (1+_dbyLine[_cl])*_font_size < 0 || _idx == av_viz.filter.sites-1) {
			_current_y = _page_h-2*_font_size;
			fprintf (av_viz.output, "showpage\n");
		}
	}
}


