angular.module('123456789', ['clotho.core', 'externalModule']);

angular.module('123456789')
.config(function() {
	console.log('widget module config block');
})
.run(function(Clotho, $timeout) {
	console.log('widget module run block');

	$timeout(function() {
		console.log('widget module timeout fn - effective callback w/o args')
	});
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
