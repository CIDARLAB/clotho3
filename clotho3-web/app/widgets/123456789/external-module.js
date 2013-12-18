angular.module('externalModule', [])
.run(function() {
	console.log('external module run block');
})
.directive('external', function () {
	return {
		restrict: 'E',
		link: function (scope, element, attrs) {
			element.html('<div class="alert alert-info">this is the external directive</div>');
		}
	};
});

