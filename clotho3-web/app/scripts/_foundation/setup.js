//adds Clotho $rootScope
angular.module('clotho.setup', [])
	.run(function ($rootScope, Clotho) {

		//extend scope with Clotho API so don't need to do this in each controller.
		$rootScope.Clotho = Clotho;

	});