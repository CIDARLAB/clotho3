
angular.module('clotho.tokenizer')

	.constant('ClothoReferenceDelimiter', {
		symbol : '@',
		keycode : 64
	})
/**
 * @ngdoc service
 * @name clothoReference
 * @description
 * Given a string, will find reference delimiters (@) and convert to links
 */
	.service('ClothoReference', function (ClothoReferenceDelimiter) {

		//wont match emails
		var matchRegexp = new RegExp('(^|[^a-zA-Z0-9-_\\.])'+ ClothoReferenceDelimiter.symbol +'([A-Za-z0-9-_\\.]+)', 'gi');

		//todo update REST url
		var urlRootREST = 'http://www.clothocad.org/data/';

		//given a function to use for wrapping a match, passes:
		// match (not include delimiter)
		// delimiter
		// index
		// textSoFar
		function convertText (wrapFn, text) {
			if (angular.isUndefined(text) || !text.length) {
				return '';
			}
			return text.replace(matchRegexp, function(match, m1, m2, index, textSoFar){
				return m1 + wrapFn(m2, ClothoReferenceDelimiter.symbol, index, textSoFar);
			});
		}

		this.convert = convertText;

		this.convertHtml = angular.bind(null, convertText, function clothoReferenceParseHtml(match) {
			return '<a sharable-popup sharable-popup-id="'+match+'" sharable-popup-trigger="mouseenter">' + match + '</a>';
		});
		
		this.convertMarkdown = angular.bind(null, convertText, function clothoReferenceParseMarkdown (match) {
			return '[' + match + '](' + urlRootREST + match + ')';
		});

		this.convertWiki = angular.bind(null, convertText, function clothoReferenceParseWiki (match) {
			return '[' + urlRootREST + match + ' ' + match + ']';
		});

	})
/**
 * @ngdoc directive
 * @name clothoReferenceParser
 * @description
 * Will go through text model and replace references with HTML (default) or Wiki, then COMPILE and set as text contents
 *
 * For now, only replaces text and not value (i.e. don't use on textarea)
 *
 * @attr clothoReferenceParser Model string to parse
 * @attr clothoReferenceType Conversion type e.g. `clotho-reference-type="wiki"`. Default is html
 *
 * @example

 <div ng-init="sample = '@hey you baby@girl whats up @clotho.developer.maxbates with ya?'">
 <div clotho-reference-parser="sample"></div>
 </div>
 *
 */
	.directive('clothoReferenceParser', function (ClothoReference, $compile, $timeout, $parse) {
		return function (scope, element, attrs) {
			scope.$watch(attrs.clothoReferenceParser, function (newval) {
				if (attrs.clothoReferenceType === 'markdown') {
					element.text(ClothoReference.convertMarkdown(newval));
				} else if (attrs.clothoReferenceType === 'wiki') {
					element.text(ClothoReference.convertWiki(newval));
				} else {
					element.html(ClothoReference.convertHtml(newval));
					console.log(element.contents());
					$compile(element.contents())(scope);
				}
			});
		};
	})
/**
 * @ngdoc filter
 *
 * @description simple filter to parse wiki or html. note that html isn't compiled. can't pass many options either
 *
 * @example

<div ng-init="sample = '@hey you baby@girl whats up @clotho.developer.maxbates with ya?'">
  <p ng-bind="sample | clothoReferenceParse:'markdown'"></p>
  <div ui-markdown="sample | clothoReferenceParse:'markdown'"></div>
</div>
 */
	.filter('clothoReferenceParse', function (ClothoReference) {
		return function (input, type, fn) {
			if (angular.isFunction(fn)) {
				return fn(input);
			}

			if (angular.isUndefined(type)) {
				return input;
			}

			switch (type) {
				case 'html' :
				{
					return ClothoReference.convertHtml(input);
				}
				case 'wiki' : {
					return ClothoReference.convertWiki(input);
				}
				case 'markdown' : {
					return ClothoReference.convertMarkdown(input);
				}
			}
		};
	});