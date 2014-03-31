/**
 @description Converts markdown into HTML
 @note relies on Showdown being present: https://github.com/coreyti/showdown

 @example
 <ui-markdown>
 #Markdown directive
 *It works!*
 </ui-markdown>

 You can also bind the markdown input to a scope variable:

 <div ui-markdown="markdown">
 </div>

 */

angular.module('clotho.interface').directive('uiMarkdown', function () {
	var converter;
	var showdownPromise = $clotho.extensions.mixin('bower_components/showdown/src/showdown.js').then(function() {
		converter = new Showdown.converter();
	});
	
	return {
		restrict: 'EA',
		link: function (scope, element, attrs) {
			showdownPromise.then(function() {
				if (attrs.uiMarkdown) {
					scope.$watch(attrs.uiMarkdown, function (newVal) {
						var html = newVal ? converter.makeHtml(newVal) : '';
						element.html(html);
					});
				} else {
					var html = converter.makeHtml(element.text());
					element.html(html);
				}
			});
		}
	};
});