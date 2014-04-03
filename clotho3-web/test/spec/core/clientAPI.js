describe('clotho.core ClientAPI', function() {

	var ClientAPI;

	// excuted before each "it()" is run.
	beforeEach(function() {
		// load the module
		module('clotho.core');

		// inject your factory for testing
		inject(function (_ClientAPI_) {
			ClientAPI = _ClientAPI_;
		});
	});

	it('should be defined', function () {
		expect(angular.isObject(ClientAPI)).toBe(true);
	});

});