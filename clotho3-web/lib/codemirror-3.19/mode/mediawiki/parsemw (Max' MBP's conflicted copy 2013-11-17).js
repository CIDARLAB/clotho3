// EXAMPLE: http://mwcodemirror.sourceforge.net/
//
// Syntax highlighting MediaWiki parser for the CodeMirror framework (http://marijn.haverbeke.nl/codemirror/)
//
// Version 0.1, October 6th, 2009
//
// Writen by Wikimedia Commons user JovanCormac (http://commons.wikimedia.org/wiki/User:JovanCormac)
//
//
// This program is free software; you can redistribute it and/or modify it under the terms of the 
// GNU General Public License as published by the Free Software Foundation.
// 
// Read the full license at http://www.opensource.org/licenses/gpl-3.0.html
//

var MWParser = Editor.Parser = (function() {
	var tokenizeMW = (function() {
		function normal(source, setState) {
			var ch = source.next();
			
			if (ch == "<" && source.lookAhead("!--", true)) {
				// Comment
				setState(inComment);
				return null;
			}
			else if (ch == "'") {
				if (source.lookAhead("''''", true)) {
					// Bold & italic text
					setState(inBoldItalic);
					return null;
				}
				else if (source.lookAhead("''", true)) {
					// Bold text
					setState(inBold);
					return null;
				}
				else if (source.lookAhead("'", true)) {
					// Italic text
					setState(inItalic);
					return null;
				}
				else {
					// Normal wikitext
					source.nextWhileMatches(/[^\n\[{<']/);
					return "mw-text";
				}
			}
			else if (ch == "[") {
				if (source.lookAhead("[File:", true, false, true)) {
					// File
					setState(inFile);
					return null;
				}
				else if (source.lookAhead("[", true)) {
					// Interwiki link
					setState(inLink);
					return null;
				}
				else {
					// Weblink
					setState(inWeblink);
					return null;
				}
			}
			else if (ch == "{") {
				if (source.lookAhead("|", true)) {
					// Table
					setState(inTable);
					return null;
				}
				else {
					// Normal wikitext
					source.nextWhileMatches(/[^\n\[{<']/);
					return "mw-text";
				}
			}
			else {
				// Normal wikitext
				source.nextWhileMatches(/[^\n\[{<']/);
				return "mw-text";
			}
		}
		
		function inComment(source, setState) {
			while (!source.endOfLine()) {
				var ch = source.next();
				if (ch == "-" && source.lookAhead("->", true)) {
					setState(normal);
					break;
				}
			}
			return "mw-comment";
		}
		
		function inBoldItalic(source, setState) {
			while (!source.endOfLine()) {
				var ch = source.next();
				if (ch == "'" && source.lookAhead("''''", true)) {
					setState(normal);
					break;
				}
			}
			return "mw-bolditalic";
		}
		
		function inBold(source, setState) {
			while (!source.endOfLine()) {
				var ch = source.next();
				if (ch == "'" && source.lookAhead("''", true)) {
					setState(normal);
					break;
				}
			}
			return "mw-bold";
		}
		
		function inItalic(source, setState) {
			while (!source.endOfLine()) {
				var ch = source.next();
				if (ch == "'" && source.lookAhead("'", true)) {
					setState(normal);
					break;
				}
			}
			return "mw-italic";
		}
		
		function inFile(source, setState) {
			var closed = false;
			
			while (!source.endOfLine()) {
				var ch = source.next();
				if (ch == "]" && source.lookAhead("]", true)) {
					closed = true;
					break;
				}
			}
			
			setState(normal);
			if (closed) return "mw-file"; else return "mw-syntaxerror";
		}
		
		function inLink(source, setState) {
			var closed = false;
			
			while (!source.endOfLine()) {
				var ch = source.next();
				if (ch == "]" && source.lookAhead("]", true)) {
					closed = true;
					break;
				}
			}
			
			setState(normal);
			if (closed) return "mw-link"; else return "mw-syntaxerror";
		}
		
		function inWeblink(source, setState) {
			var closed = false;
			
			while (!source.endOfLine()) {
				var ch = source.next();
				if (ch == "]") {
					closed = true;
					break;
				}
			}
			
			setState(normal);
			if (closed) return "mw-weblink"; else return "mw-syntaxerror";
		}
		
		function inTable(source, setState) {
			while (!source.endOfLine()) {
				var ch = source.next();
				if (ch == "|" && source.lookAhead("}", true)) {
					setState(normal);
					break;
				}
			}
			return "mw-table";
		}
		
		return function(source, startState) {
			return tokenizer(source, startState || normal);
		};
	})();
	
	function parseMW(source, space) {
		function indentTo(n) {return function() {return n;}}
		
		var tokens = tokenizeMW(source);		
		var space = space || 0;
		
		var iter = {
			next: function() {
		        var token = tokens.next(), style = token.style, content = token.content;
				if (content == "\n") {
					token.indentation = indentTo(space);
				}
				return token;
			},
			copy: function() {
				var _tokenState = tokens.state;
				return function(source) {
					tokens = tokenizeMW(source, _tokenState);
					return iter;
				};
			}
		};
		return iter;
	}
	return {make: parseMW};
})();
