describe('clotho.core PubSub', function() {

	var PubSub;

	// excuted before each "it()" is run.
	beforeEach(function() {
		// load the module
		module('clotho.core');

		// inject your factory for testing
		inject(function (_PubSub_) {
			PubSub = _PubSub_;
		});
	});

	it('should be defined and have functions: trigger, on, off, once, destroy, clear', function() {
		expect(angular.isObject(PubSub)).toBe(true);
		expect(angular.isFunction(PubSub.trigger)).toBe(true);
		expect(angular.isFunction(PubSub.on)).toBe(true);
		expect(angular.isFunction(PubSub.once)).toBe(true);
		expect(angular.isFunction(PubSub.off)).toBe(true);
		expect(angular.isFunction(PubSub.destroy)).toBe(true);
		expect(angular.isFunction(PubSub.clear)).toBe(true);
	});


	it('', function () {});
	it('', function () {});
	it('', function () {});
	it('', function () {});
	it('', function () {});
	it('', function () {});

});