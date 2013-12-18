angular.module('123456789', ['clotho.core', 'externalModule']);

angular.module('123456789').run(function(Clotho) {
	console.log('logging clotho shows dependencies can be shared from the main app');
	console.log(Clotho);
	console.log('this module also depends on externalModule, which should log something in the console');
});

angular.module('123456789').directive('special', function () {

	return {
		restrict: 'E',
		template: '<div class="well">' +
			'<p>this is the directive from the new module</p>' +
			'<p>This value is from parent (outside app - won\'t interpolate): {{ someValue }}</p>' +
			'<p>This interpolated value is from widget: {{ specialValue }}</p>' +
			'<p>Below is a directive that relies on externalModule</p>' +
			'<external><external>' +
		'</div>',
		link: function (scope, element, attrs) {
			scope.specialValue = 'this is from the widget';
		}
	};
});
