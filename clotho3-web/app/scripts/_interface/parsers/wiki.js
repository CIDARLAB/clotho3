'use strict';

/*
use as element to parse interally
use as attribute with model to parse dynamically based on variable (e.g. ngModel)
conversion based on http://en.wikipedia.org/wiki/Wikipedia:Cheatsheet,
however, interpolation of {{ }}  templating is left to angular
 */
angular.module('clotho.interface')
  .directive('wiki', function () {

		function wiki2html(s) {
			
			console.log(s);
			
			if (!s) return s;

			// lists need to be done using a function to allow for recusive calls
			function list(str) {
				return str.replace(/(?:(?:(?:^|\n)[\*#].*)+)/g, function (m) {  // (?=[\*#])
					var type = m.match(/(^|\n)#/) ? 'OL' : 'UL';
					// strip first layer of list
					m = m.replace(/(^|\n)[\*#][ ]{0,1}/g, "$1");
					m = list(m);
					return '<' + type + '><li>' + m.replace(/^\n/, '').split(/\n/).join('</li><li>') + '</li></' + type + '>';
				});
			}

			return list(s

				/* BLOCK ELEMENTS */
				.replace(/(?:^|\n+)([^# =\*<].+)(?:\n+|$)/gm, function (m, l) {
					if (l.match(/^\^+$/)) return l;
					return "\n<p>" + l + "</p>\n";
				})

				.replace(/(?:^|\n)[ ]{2}(.*)+/g, function (m, l) { // blockquotes
					if (l.match(/^\s+$/)) return m;
					return '<blockquote>' + l + '</pre>';
				})

				.replace(/((?:^|\n)[ ]+.*)+/g, function (m) { // code
					if (m.match(/^\s+$/)) return m;
					return '<pre>' + m.replace(/(^|\n)[ ]+/g, "$1") + '</pre>';
				})

				.replace(/(?:^|\n)([=]+)(.*)\1/g, function (m, l, t) { // headings
					return '<h' + l.length + '>' + t + '</h' + l.length + '>';
				})

				/* INLINE ELEMENTS */
				.replace(/'''(.*?)'''/g, function (m, l) { // bold
					return '<strong>' + l + '</strong>';
				})

				.replace(/''(.*?)''/g, function (m, l) { // italic
					return '<em>' + l + '</em>';
				})

				.replace(/[^\[](http[^\[\s]*)/g, function (m, l) { // normal link
					return '<a href="' + l + '">' + l + '</a>';
				})

				.replace(/[\[](http.*)[!\]]/g, function (m, l) { // external link
					var p = l.replace(/[\[\]]/g, '').split(/ /);
					var link = p.shift();
					return '<a href="' + link + '">' + (p.length ? p.join(' ') : link) + '</a>';
				})

				.replace(/\[\[(.*?)\]\]/g, function (m, l) { // internal link or image
					var p = l.split(/\|/);
					var link = p.shift();

					if (link.match(/^Image:(.*)/)) {
						// no support for images - since it looks up the source from the wiki db :-(
						return m;
					} else {
						return '<a href="' + link + '">' + (p.length ? p.join('|') : link) + '</a>';
					}
				})
			);
		}

    return {
      template: '<div></div>',
	    replace : true,
      restrict: 'EA',
	    scope: {},
      link: function postLink(scope, element, attrs) {
	      if (attrs.wiki) {
		      scope.$watch(attrs.wiki, function (newval, oldval) {
			      if (oldval != newval)
				      element.html(wiki2html(newval));
		      })
	      } else {
		      element.html(wiki2html(element.text()));
	      }

      }
    };
  });
