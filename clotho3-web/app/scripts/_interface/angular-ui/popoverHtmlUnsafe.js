// Angular UI Bootstrap provide some excellent directives, but the popover didn't allow for HTML content
// The popoverHtmlUnsafe and popoverHtmlUnsafePopup implement this on top of the AngularUI Bootstrap's $tooltip service
angular.module('clotho.interface').directive( 'popoverHtmlUnsafePopup', function ($templateCache) {

	$templateCache.put("template/popover/popover-html-unsafe-popup.html",
			"<div class=\"popover {{placement}}\" ng-class=\"{ in: isOpen(), fade: animation() }\">\n" +
			"  <div class=\"arrow\"></div>\n" +
			"\n" +
			"  <div class=\"popover-inner\">\n" +
			"      <h3 class=\"popover-title\" ng-bind=\"title\" ng-show=\"title\"></h3>\n" +
			"      <div class=\"popover-content\" bind-html-unsafe=\"content\"></div>\n" +
			"  </div>\n" +
			"</div>\n" +
			"");

	return {
		restrict: 'EA',
		replace: true,
		scope: { title: '@', content: '@', placement: '@', animation: '&', isOpen: '&' },
		templateUrl: 'template/popover/popover-html-unsafe-popup.html'
	};
})

	.directive( 'popoverHtmlUnsafe', function ( $compile, $timeout, $parse, $window, $tooltip ) {
		return $tooltip( 'popoverHtmlUnsafe', 'popover', 'click' );
	});