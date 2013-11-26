angular.module('clotho.interface').directive('closeable', function($compile) {
	return {
		restrict : 'A',
		priority : 1100,
		link: function compile(scope, element, attrs) {
			scope.removeElement = function() {
				element.remove();
			};
			element.prepend($compile('<a class="close" style="position: absolute; top: 12px; right: 15px;" ng-click="removeElement()">&times;</a>')(scope));
		}
	}
});