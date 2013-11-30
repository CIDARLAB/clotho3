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

Application.Interface.directive('uiMarkdown', [function () {
		var converter = new Showdown.converter();
		return {
			restrict: 'EA',
			link: function (scope, element, attrs) {
				if (attrs.uiMarkdown) {
					scope.$watch(attrs.uiMarkdown, function (newVal) {
						var html = newVal ? converter.makeHtml(newVal) : '';
						element.html(html);
					});
				} else {
					var html = converter.makeHtml(element.text());
					element.html(html);
				}
			}
		};
	}]);